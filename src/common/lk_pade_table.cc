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

#include "common/lk_pade_table.h"

#include <stdint.h>
#include <stdio.h>

#include <cassert>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/unordered_map.hpp>

#include "common/conf.h"
#include "common/rho_finder.h"
#include "common/version_number.h"
#include "find_confs/find_confs.h"

LkTable::LkTable()
    : version_bit_string_(),
      version_number_(),
      num_confs_(),
      theta_(),
      interpolate_() { }

PadeTable::PadeTable()
    : version_bit_string_(),
      version_number_(),
      num_confs_(),
      theta_(),
      num_coeffs_(),
      defect_threshold_() { }

PadeComputation::PadeComputation(LkTable *lk_table,
                                 PadeTable const &pade_table,
                                 ConfPadeMap::const_iterator start_iter,
                                 ConfPadeMap::const_iterator end_iter)
    : lk_table_(lk_table),
      pade_table_(pade_table),
      conf_lk_map_(lk_table->conf_lk_map_),
      conf_pade_map_(pade_table.conf_pade_map_),
      start_iter_(start_iter),
      end_iter_(end_iter) { }

double PadeComputation::ComputePadeLk(PadeCoefficients const &pade_coefficients,
                                      double pade_defect_threshold,
                                      double rho) {
  std::vector<double> const &coeffs = pade_coefficients.coeffs_;
  std::vector<std::vector<double> > const &roots = pade_coefficients.roots_;

  uint64_t num_coeffs_to_remove = 0;
  for (size_t coeff_set_id = 0; coeff_set_id < roots.size(); ++coeff_set_id) {
    bool root_nearby = false;
    std::vector<double> const &roots_for_coeff = roots[coeff_set_id];
    for (size_t root_id = 0; root_id < roots_for_coeff.size(); ++root_id) {
      double root = roots_for_coeff[root_id];
      if (std::abs(rho - root) < pade_defect_threshold) {
        root_nearby = true;
        break;
      }
    }

    // coeff_set_id is equal to the number of coefficients to remove.
    num_coeffs_to_remove = coeff_set_id;

    // No root nearby, so no need to remove any more coefficients.
    if (!root_nearby) {
      break;
    }

    // If coeff_set_id is final coefficient _with_ roots, and there is still
    // a root nearby, then increment num_coeffs_to_remove.
    // This will then be sufficient to avoid all roots.
    if (coeff_set_id == roots.size() - 1 && root_nearby) {
      ++num_coeffs_to_remove;
    }
  }

  assert(num_coeffs_to_remove < coeffs.size());

  uint64_t num_coeffs_used = coeffs.size() - num_coeffs_to_remove;

  assert(num_coeffs_used > 0);

  // Use Pade coefficients to compute the likelihood for conf.
  double denominator = 1.0;
  for (uint64_t coeff_id = 0; coeff_id < num_coeffs_used - 1; ++coeff_id) {
    denominator =
      1.0 + coeffs[num_coeffs_used - 1 - coeff_id] * (1.0/rho) / denominator;
  }
  double lk = coeffs[0] / denominator;

  assert(lk <= 1.0);

  return lk;
}

void PadeComputation::operator()() {
  for (ConfPadeMap::const_iterator conf_iter = start_iter_;
    conf_iter != end_iter_;
    ++conf_iter) {
    std::vector<double> &lk_column =
      conf_lk_map_.at(conf_iter->first).second;
    for (size_t rho_id = lk_table_->num_lk_rhos_;
      rho_id < lk_table_->rho_finder_.rho_list_.size();
      ++rho_id) {
      double lk = ComputePadeLk(conf_iter->second.second,
                                pade_table_.defect_threshold_,
                                lk_table_->rho_finder_.rho_list_[rho_id]);
      // If lk is not positive, use truncated likelihood.
      if (lk <= 0.0) {
        lk = lk_column.at(lk_table_->num_lk_rhos_ - 1);
      }
      lk_column.at(rho_id) = lk;
    }
  }
}

LkTable LoadLikelihoodAndPade(std::string const &lk_file,
                              std::string const &pade_file,
                              double const pade_resolution,
                              double const pade_max,
                              uint32_t const window_size,
                              std::vector<std::string> const &snp_seqs,
                              uint32_t const num_threads) {
  printf("Loading likelihood and Pade files.\n");

  if (snp_seqs.size() == 0) {
    fprintf(stderr, "There must be at least one sequence.\n");
    std::exit(1);
  };

  if (window_size < 2) {
    fprintf(stderr, "Window size must be at leaset 2.\n");
    std::exit(1);
  }

  // Find confs in data.
  assert(snp_seqs.size() >= 2);
  size_t len_snp_seq = snp_seqs.front().size();
  assert(len_snp_seq >= 2);

  std::set<Conf> confs_set = FindConfsMT(num_threads, window_size, snp_seqs);
  std::vector<Conf> confs(confs_set.begin(), confs_set.end());

  // Read data from file.
  PadeTable pade_table;
  LkTable lk_table;

  ConfLkMap &conf_lk_map = lk_table.conf_lk_map_;
  ConfPadeMap &conf_pade_map = pade_table.conf_pade_map_;

  int fread_ret;

  FILE *fp_lk = fopen(lk_file.c_str(), "rb");
  if (fp_lk == NULL) {
    fprintf(stderr,
            "Unable to open likelihood file for reading: %s.\n",
            lk_file.c_str());
    std::exit(1);
  }

  fread_ret = fread(reinterpret_cast<void *>(&lk_table.version_bit_string_),
                   sizeof(lk_table.version_bit_string_), 1, fp_lk);
  assert(fread_ret == 1);

  lk_table.version_bit_string_ = lk_table.version_bit_string_ - TABLE_GEN_SALT;
  if (lk_table.version_bit_string_ != TABLE_GEN_VERSION_OUTPUT) {
    fprintf(stderr, "The likelihood file has the wrong version number.\n");
    std::exit(1);
  }

  fread_ret = fread(reinterpret_cast<void *>(&lk_table.num_confs_),
                    sizeof(lk_table.num_confs_), 1, fp_lk);
  assert(fread_ret == 1);
  fread_ret = fread(reinterpret_cast<void *>(&lk_table.theta_),
                    sizeof(lk_table.theta_), 1, fp_lk);
  assert(fread_ret == 1);
  printf("Theta = %f\n", lk_table.theta_);
  fread_ret = fread(reinterpret_cast<void *>(&lk_table.interpolate_),
                    sizeof(lk_table.interpolate_), 1, fp_lk);
  assert(fread_ret == 1);
  if(lk_table.interpolate_){
     printf("We will use linear interpolation for this table.\n");
  } else {
     printf("We will not use linear interpolation for this table.\n");
  }

  FILE *fp_pade = NULL;

  if (!pade_file.empty()) {
    fp_pade = fopen(pade_file.c_str(), "rb");
    if (fp_pade == NULL) {
      fprintf(stderr,
              "Unable to open Pade file for reading: %s.\n",
              pade_file.c_str());
      std::exit(1);
    }
    fread_ret = fread(reinterpret_cast<void *>(&pade_table.version_bit_string_),
                      sizeof(pade_table.version_bit_string_), 1, fp_pade);
    assert(fread_ret == 1);
    pade_table.version_number_ = pade_table.version_bit_string_ - PADE_SALT;

    if (pade_table.version_number_ != PADE_VERSION_OUTPUT) {
      fprintf(stderr, "The Pade file has the wrong version number.\n");
      std::exit(1);
    }

    fread_ret = fread(reinterpret_cast<void *>(&pade_table.num_confs_),
                      sizeof(pade_table.num_confs_), 1, fp_pade);
    assert(fread_ret == 1);
    fread_ret = fread(reinterpret_cast<void *>(&pade_table.theta_),
                      sizeof(pade_table.theta_), 1, fp_pade);
    assert(fread_ret == 1);
    fread_ret = fread(reinterpret_cast<void *>(&pade_table.num_coeffs_),
                      sizeof(pade_table.num_coeffs_), 1, fp_pade);
    assert(fread_ret == 1);
    fread_ret = fread(reinterpret_cast<void *>(&pade_table.defect_threshold_),
                      sizeof(pade_table.defect_threshold_), 1, fp_pade);
    assert(fread_ret == 1);

    if (pade_table.theta_ != lk_table.theta_) {
      fprintf(stderr,
              "The likelihood file has a different theta from "
              "the pade file.\n.");
      std::exit(1);
    }
  }

  uint64_t num_rho_segments;
  fread_ret = fread(reinterpret_cast<void *>(&num_rho_segments),
                    sizeof(num_rho_segments), 1, fp_lk);
  assert(fread_ret == 1);
  {
    std::vector<double> rho_range;
    for (size_t i = 0; i < num_rho_segments; ++i) {
      double start, delta;
      fread_ret = fread(reinterpret_cast<void *>(&start),
                        sizeof(start), 1, fp_lk);
      assert(fread_ret == 1);
      fread_ret = fread(reinterpret_cast<void *>(&delta),
                        sizeof(delta), 1, fp_lk);
      assert(fread_ret == 1);
      rho_range.push_back(start);
      rho_range.push_back(delta);
    }
    double end;
    fread_ret = fread(reinterpret_cast<void *>(&end), sizeof(end), 1, fp_lk);
    assert(fread_ret == 1);
    rho_range.push_back(end);

    // Compute number of rhos in the likelihood table.
    lk_table.num_lk_rhos_ =
      RhoFinder(ParseRhoRange(rho_range)).rho_list_.size();

    // Add rho values for Pade approximants.
    if (!pade_file.empty()) {
      if (pade_max > end) {
        rho_range.push_back(pade_resolution);
        rho_range.push_back(pade_max);
      }
    }

    lk_table.rho_finder_ = RhoFinder(ParseRhoRange(rho_range));
  }

  // Allocate enough space for all rho values, even those
  // reserved for Pade approximants.
  for (std::vector<Conf>::const_iterator citer = confs.begin();
       citer != confs.end();
       ++citer) {
    conf_lk_map[*citer] = std::make_pair(
        -1, std::vector<double>(lk_table.rho_finder_.rho_list_.size()));

    // Allocate space for pade coefficients.
    if (!pade_file.empty()) {
      conf_pade_map[*citer] =
        std::make_pair(-1, PadeCoefficients(pade_table.num_coeffs_));
    }
  }

  // Read in confs from likelihood file.
  {
    std::vector<Conf>::const_iterator conf_iter = confs.begin();
    for (uint64_t conf_id = 0; conf_id < lk_table.num_confs_; ++conf_id) {
      Conf::BinaryRep conf_binary_rep;
      fread_ret = fread(reinterpret_cast<void *>(&conf_binary_rep.rep_),
                        sizeof(conf_binary_rep.rep_), 1, fp_lk);
      assert(fread_ret == 1);
      Conf conf(conf_binary_rep);
      assert(conf.ComputeDegree() >= 1);

      if (*conf_iter == conf) {
        conf_lk_map[*conf_iter].first = conf_id;

        ++conf_iter;
        if (conf_iter == confs.end()) {
          fseek(fp_lk,
                (lk_table.num_confs_ - 1 - conf_id)
                  * sizeof(conf_binary_rep.rep_),
                SEEK_CUR);
          break;
        }
      }
    }
    if (conf_iter != confs.end()) {
      fprintf(stderr,
              "Configuration in data not found in likelihood file "
              "while reading in configurations.\n"
              "Configuration: %s.\n", conf_iter->GetString().c_str());
      std::exit(1);
    }
  }

  // Read in confs from pade file.
  if (!pade_file.empty()) {
    std::vector<Conf>::const_iterator conf_iter = confs.begin();
    for (uint64_t conf_id = 0; conf_id < pade_table.num_confs_; ++conf_id) {
      Conf::BinaryRep conf_binary_rep;
      fread_ret = fread(reinterpret_cast<void *>(&conf_binary_rep.rep_),
                        sizeof(conf_binary_rep.rep_), 1, fp_pade);
      assert(fread_ret == 1);
      Conf conf(conf_binary_rep);
      assert(conf.ComputeDegree() >= 1);

      if (*conf_iter == conf) {
        conf_pade_map[*conf_iter].first = conf_id;
        ++conf_iter;
        if (conf_iter == confs.end()) {
          fseek(fp_pade,
                (pade_table.num_confs_ - 1 - conf_id)
                  * sizeof(conf_binary_rep.rep_),
                SEEK_CUR);
          break;
        }
      }
    }

    if (conf_iter != confs.end()) {
      fprintf(stderr,
              "Configuration in data not found in Pade file "
              "while reading in configurations\n."
              "Configuration: %s.\n", conf_iter->GetString().c_str());
      std::exit(1);
    }
  }

  printf("Number of rho values in likelihood table: %d\n",
         static_cast<int>(lk_table.num_lk_rhos_));
  printf("Reading likelihoods: ");
  fflush(stdout);

  // Read in lks.
  for (int rho_id = 0;
       rho_id < static_cast<int>(lk_table.num_lk_rhos_);
       ++rho_id) {
    printf("%d ", rho_id);
    fflush(stdout);

    std::vector<Conf>::const_iterator conf_iter = confs.begin();
    for (uint64_t conf_id = 0; conf_id < lk_table.num_confs_; ++conf_id) {
      double prob;
      fread_ret = fread(reinterpret_cast<void *>(&prob),
                        sizeof(prob), 1, fp_lk);
      assert(fread_ret == 1);
      if (feof(fp_lk) || ferror(fp_lk)) {
        fprintf(stderr,
                "Error reading likelihood table. "
                "Possibly not enough likelihoods in the likelihood table.\n");
        std::exit(1);
      }

      if (prob < 0.0 || prob > 1.0) {
        fprintf(stderr,
               "Likelihoods must be non-negative and less than 1.0.\n");
        std::exit(1);
      }

      if (conf_lk_map[*conf_iter].first == conf_id) {
        conf_lk_map[*conf_iter].second[rho_id] = prob;
        ++conf_iter;
        if (conf_iter == confs.end()) {
          fseek(fp_lk,
                (lk_table.num_confs_ - 1 - conf_id) * sizeof(prob),
                SEEK_CUR);
          break;
        }
      }
    }

    if (conf_iter != confs.end()) {
      fprintf(stderr,
              "Configuration in data not found in likelihood table "
              "while reading in likelihoods.\n"
              "Configuration: %s.\n", conf_iter->GetString().c_str());
      std::exit(1);
    }
  }

  printf("\n");

  assert(!feof(fp_lk) && !ferror(fp_lk));

  fgetc(fp_lk);

  if (!feof(fp_lk)) {
    fprintf(stderr, "Too many likelihoods in likelihood table or "
                    "likelihood file is corrupt.\n");
    std::exit(1);
  }

  assert(fp_lk != NULL);
  fclose(fp_lk);

  if (!pade_file.empty()) {
    std::vector<Conf>::const_iterator conf_iter = confs.begin();

    // Read in Pade coefficients.
    for (uint64_t conf_id = 0; conf_id < pade_table.num_confs_; ++conf_id) {
      bool relevant_conf;
      if (conf_iter != confs.end()) {
        relevant_conf = conf_pade_map[*conf_iter].first == conf_id;
      } else {
        relevant_conf = false;
      }

      for (uint64_t coeff_id = 0;
           coeff_id < pade_table.num_coeffs_;
           ++coeff_id) {
        double pade_coeff;
        fread_ret = fread(reinterpret_cast<void *>(&pade_coeff),
                          sizeof(pade_coeff), 1, fp_pade);
        assert(fread_ret == 1);
        if (feof(fp_pade) || ferror(fp_pade)) {
          fprintf(stderr,
                  "Error reading Pade file. "
                  "Possibly not enough Pade coefficients in the Pade file.\n");
          std::exit(1);
        }
        if (relevant_conf) {
          conf_pade_map[*conf_iter].second.coeffs_[coeff_id] = pade_coeff;
        }
      }

      // Read in Pade roots.
      // Read in number of sets of coefficients with roots.
      uint8_t num_coeff_sets_from_file;
      fread_ret = fread(reinterpret_cast<void *>(&num_coeff_sets_from_file),
                        sizeof(num_coeff_sets_from_file), 1, fp_pade);
      assert(fread_ret == 1);
      if (feof(fp_pade) || ferror(fp_pade)) {
        fprintf(stderr,
                "Error reading Pade file. "
                "Possibly not enough Pade coefficients in the Pade file.\n");
        std::exit(1);
      }
      uint32_t num_coeff_sets = static_cast<uint32_t>(num_coeff_sets_from_file);

      if (relevant_conf) {
        conf_pade_map[*conf_iter].second.roots_.resize(num_coeff_sets);
      }

      for (uint32_t coeff_set_id = 0;
           coeff_set_id < num_coeff_sets;
           ++coeff_set_id) {
        uint8_t num_roots_from_file;
        fread_ret = fread(reinterpret_cast<void *>(&num_roots_from_file),
                          sizeof(num_roots_from_file), 1, fp_pade);
        assert(fread_ret == 1);
        uint32_t num_roots = static_cast<uint32_t>(num_roots_from_file);

        assert(num_roots < pade_table.num_coeffs_);

        if (relevant_conf) {
          conf_pade_map[*conf_iter].second.roots_[coeff_set_id]
                                             .resize(num_roots);
        }

        for (uint32_t root_id = 0; root_id < num_roots; ++root_id) {
          double root;
          fread_ret = fread(reinterpret_cast<void *>(&root),
                            sizeof(root), 1, fp_pade);
          assert(fread_ret == 1);
          if (relevant_conf) {
            conf_pade_map[*conf_iter].second
              .roots_[coeff_set_id][root_id] = root;
          }
        }
      }

      if (relevant_conf) {
        ++conf_iter;
      }
    }

    assert(!feof(fp_pade) && !ferror(fp_pade));

    fgetc(fp_pade);

    if (!feof(fp_pade)) {
      fprintf(stderr, "Too many values in Pade file or "
                      "Pade file is corrupt.\n");
      std::exit(1);
    }

    if (conf_iter != confs.end()) {
      fprintf(stderr,
              "Configuration in data not found in Pade file "
              "while reading in Pade coefficients.\n"
              "Configuration: %s.\n", conf_iter->GetString().c_str());
    }

    assert(fp_pade != NULL);
    fclose(fp_pade);

    {
      printf("Precomputing Pade likelihoods.\n");

      uint64_t partition_length = conf_pade_map.size()/num_threads;

      std::vector<std::pair<ConfPadeMap::const_iterator,
                            ConfPadeMap::const_iterator> > partitions;
      ConfPadeMap::const_iterator c_iter = conf_pade_map.begin();
      assert(num_threads > 0);
      for (uint32_t partition_id = 0;
           partition_id < num_threads - 1;
           ++partition_id) {
        ConfPadeMap::const_iterator start_iter = c_iter;
        std::advance(c_iter, partition_length);
        ConfPadeMap::const_iterator end_iter = c_iter;
        partitions.push_back(std::make_pair(start_iter, end_iter));
      }

      partitions.push_back(std::make_pair(c_iter, conf_pade_map.end()));

      assert(partitions.size() == num_threads);

      boost::thread_group tgroup;
      for (uint32_t partition_id = 0;
           partition_id < partitions.size();
           ++partition_id) {
        ConfPadeMap::const_iterator start_iter, end_iter;
        boost::tie(start_iter, end_iter) = partitions[partition_id];
        tgroup.create_thread(
            PadeComputation(&lk_table, pade_table, start_iter, end_iter));
      }
      tgroup.join_all();
    }
  }

  return lk_table;
}
