// Copyright (C) 2012  Andrew H. Chan, Paul A. Jenkins, Yun S. Song
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// <http://www.gnu.org/licenses/>
//
// 

#include "common/mut_mat_prior.h"

#include <stdint.h>

#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

MutationMatrix NormalizedMutMatrix(MutationMatrix const &mut_mat) {
  MutationMatrix norm_mut_mat;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      norm_mut_mat[i][j] = 0.0;
    }
  }

  // Normalize.
  for (int i = 0; i < 4; ++i) {
    double norm = 0.0;
    for (int j = 0; j < 4; ++j) {
      norm += mut_mat[i][j];
    }
    for (int j = 0; j < 4; ++j) {
      norm_mut_mat[i][j] = mut_mat[i][j]/norm;
    }
  }

  return norm_mut_mat;
}

MutationMatrix LoadMutationMatrixHelper(std::string const &mut_mat_file) {
  // Load mutation matrix.
  // {ACGT} x {ACGT}; row indicates ancestor, column indicates derived.
  FILE *fp = fopen(mut_mat_file.c_str(), "r");
  if (fp == NULL) {
    fprintf(stderr,
            "Unable to open mutation matrix file: %s.\n", mut_mat_file.c_str());
    exit(1);
  }

  MutationMatrix mut_mat;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      int number_read = fscanf(fp, "%lf", &mut_mat[i][j]);
      if (feof(fp) || ferror(fp) || number_read != 1) {
        fprintf(stderr, "Mutation matrix file has too few elements.\n");
        exit(1);
      }
    }
  }

  assert(!feof(fp) && !ferror(fp));
  int number_read = fscanf(fp, "%*f");
  if (number_read != -1) {
    fprintf(stderr, "Mutation matrix file has too many elements.\n");
    exit(1);
  }

  return NormalizedMutMatrix(mut_mat);
}

MutationMatrix DefaultMutationMatrix() {
    MutationMatrix mut_mat;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i == j) {
              mut_mat[i][j]=0.0;
            } else {
              mut_mat[i][j]=1.0;
            }
        }
    }
    return NormalizedMutMatrix(mut_mat);
}

MutationMatrix LoadMutationMatrix(std::string const &mut_mat_file) {
  MutationMatrix mut_mat;
  if (mut_mat_file != "") {
    mut_mat = LoadMutationMatrixHelper(mut_mat_file);
  } else {
    mut_mat = DefaultMutationMatrix();
  }

  // Display mutation matrix.
  printf("Mutation matrix:\n");
  for (int i = 0; i < static_cast<int>(mut_mat.size()); ++i) {
    assert(mut_mat[i].size() > 0);
    printf("[");
    for (int j = 0; j < static_cast<int>(mut_mat[i].size()) - 1; ++j) {
      printf("%.4f, ", mut_mat[i][j]);
    }
    printf("%.4f]\n", mut_mat[i].back());
  }

  return mut_mat;
}

Prior LoadPriorHelper(std::string const &prior_file_name,
                std::vector<uint64_t> const &snp_pos) {
  FILE *fp = fopen(prior_file_name.c_str(), "r");
  if (fp == NULL) {
    fprintf(stderr,
            "Unable to open prior file: %s.\n", prior_file_name.c_str());
    exit(1);
  }

  Prior prior;
  std::vector<uint64_t>::const_iterator pos_iter = snp_pos.begin();

  while (pos_iter != snp_pos.end()) {
    // Parse snp position and prior.
    uint32_t snp_pos1;
    PriorSite cur_prior;
    int number_read = fscanf(fp, "%d %lf %lf %lf %lf",
                             &snp_pos1, &cur_prior[0], &cur_prior[1],
                             &cur_prior[2], &cur_prior[3]);
    if (number_read != 5) {
      fprintf(stderr,
              "Problem reading in ancestral prior entries.\n");
      exit(1);
    }

    if (cur_prior[0] < 0 ||
      cur_prior[1] < 0 ||
      cur_prior[2] < 0 ||
      cur_prior[3] < 0) {
      fprintf(stderr,
              "Prior must have non-negative entries.\n");
      exit(1);
    }

    // Normalize prior.
    double prior_norm = cur_prior[0]
                      + cur_prior[1]
                      + cur_prior[2]
                      + cur_prior[3];
    cur_prior[0] /= prior_norm;
    cur_prior[1] /= prior_norm;
    cur_prior[2] /= prior_norm;
    cur_prior[3] /= prior_norm;

    if (snp_pos1 > *pos_iter) {
      fprintf(stderr,
              "The ancestral prior file is missing the prior for a SNP "
              "in the sequence file, or the ancestral priors are not "
              "ordered in increasing genomic position.");
      exit(1);
    }

    if (snp_pos1 == *pos_iter) {
      prior.push_back(cur_prior);
      ++pos_iter;
    }
  }

  assert(prior.size() == snp_pos.size());

  return prior;
}

Prior DefaultPrior(MutationMatrix const &mut_mat,
                   std::vector<uint64_t> const &snp_pos) {
    Prior prior(snp_pos.size());

    int mut_dim = mut_mat.size();
    assert(mut_mat.size() == mut_mat[0].size());

    // data is the _transpose_ of mut_mat.
    double *data = new double[mut_dim * mut_dim];
    for (size_t i = 0; i < mut_mat.size(); ++i) {
      for (size_t j = 0; j < mut_mat.size(); ++j) {
        data[i * mut_mat.size() + j] = mut_mat[j][i];
      }
    }

    // Compute dominant eigenvector.
    gsl_matrix_view gslMatrix = gsl_matrix_view_array(data, mut_dim, mut_dim);

    gsl_vector_complex *eigen_values = gsl_vector_complex_alloc(mut_dim);
    gsl_matrix_complex *eigen_vectors =
      gsl_matrix_complex_alloc(mut_dim, mut_dim);

    gsl_eigen_nonsymmv_workspace *gsl_workspace =
        gsl_eigen_nonsymmv_alloc(mut_dim);

    gsl_eigen_nonsymmv(&gslMatrix.matrix,
                       eigen_values,
                       eigen_vectors,
                       gsl_workspace);

    gsl_eigen_nonsymmv_free(gsl_workspace);

    gsl_eigen_nonsymmv_sort(eigen_values,
                            eigen_vectors,
                            GSL_EIGEN_SORT_ABS_DESC);

    gsl_vector_complex_view dom_eigen_vector =
      gsl_matrix_complex_column(eigen_vectors, 0);

    gsl_complex dom_eigen_value = gsl_vector_complex_get(eigen_values, 0);
    assert(GSL_REAL(dom_eigen_value) > 0.99 &&
           GSL_REAL(dom_eigen_value) < 1.01);

    PriorSite stationary_dist;
    for (int i = 0; i < mut_dim; ++i) {
        gsl_complex elem = gsl_vector_complex_get(&dom_eigen_vector.vector, i);
        stationary_dist[i] = GSL_REAL(elem);
    }

    // Normalize.
    double stationary_sum = 0.0;
    for (int i = 0; i < mut_dim; ++i) {
        stationary_sum += stationary_dist[i];
    }

    for (int i = 0; i < mut_dim; ++i) {
        stationary_dist[i] /= stationary_sum;
    }

    gsl_matrix_complex_free(eigen_vectors);
    gsl_vector_complex_free(eigen_values);
    delete data;

    std::fill(prior.begin(), prior.end(), stationary_dist);

    return prior;
}

Prior LoadPrior(MutationMatrix const &mut_mat,
                std::vector<uint64_t> const &snp_pos,
                std::string const &prior_file) {
  Prior prior;
  if (!prior_file.empty()) {
    printf("Loading ancestral priors from file.\n");

    prior = LoadPriorHelper(prior_file, snp_pos);
  } else {
    printf("Using stationary distribution of mutation matrix "
           "for ancestral priors.\n");

    prior = DefaultPrior(mut_mat, snp_pos);
  }

  return prior;
}
