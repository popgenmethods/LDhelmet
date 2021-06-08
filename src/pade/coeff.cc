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

#include "pade/coeff.h"

#include <stddef.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <gsl/gsl_poly.h>

#include "common/ncr.h"
#include "pade/subtable.h"

int ComputeCombinatorialFactor(Conf const &c, Conf const &r) {
  int prod = 1;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      prod *= Ncr(c.Getr(i, j), r.Getr(i, j));
    }
  }
  return prod;
}

double ComputeQValueIncrementally(Table const &table, int cur_level,
                                  int m, Conf const &c) {
  int const c_tot = c.ab_ + c.aB_ + c.Ab_ + c.AB_;
  if (m > c_tot) {
    // Add zero because the combinatorial factor is zero in this case.
    return 0.0;
  }

  double sum = 0.0;
  int const u = cur_level - m;

  ConfGenMC conf_gen_m_c(c, m);

  // Iterate over all subconfigurations r of c such that
  // r.ab_ + r.aB_ + r.Ab_ + r.AB_ == m.
  while (!conf_gen_m_c.end()) {
    Conf r = *conf_gen_m_c;
    int const combinatorial_factor = ComputeCombinatorialFactor(c, r);
    sum += combinatorial_factor * GetTable(table, m, u, r);
    ++conf_gen_m_c;
  }

  return sum;
}

Coeffs ComputeCoeffsForConf(uint32_t num_coeffs,
                            std::vector<double> const &q_values_for_conf,
                            Conf const &conf) {
    Coeffs coeffs;

    coeffs.cont_frac = std::vector<double>(num_coeffs, 0.0);
    coeffs.poly_numerator = std::vector<std::vector<double> >(
        num_coeffs,
        std::vector<double>(num_coeffs, 0.0));
    coeffs.poly_denominator = std::vector<std::vector<double> >(
        num_coeffs,
        std::vector<double>(num_coeffs, 0.0));

    // Compute A and B.
    std::vector<std::vector<double> > A(num_coeffs,
                                        std::vector<double>(num_coeffs, 0.0));
    std::vector<std::vector<double> > B(num_coeffs,
                                        std::vector<double>(num_coeffs, 0.0));

    for (uint32_t cur_level = 0; cur_level < 2*num_coeffs; cur_level += 2) {
      A[0][cur_level/2] = q_values_for_conf[cur_level/2];
      B[0][cur_level/2] = (cur_level == 0 ? 1.0 : 0.0);
      for (uint32_t i = 2; i <= cur_level; i += 2) {
        A[i/2][(cur_level - i)/2] =
          B[i/2 - 1][(cur_level - i)/2 + 1]
          - A[i/2 - 1][(cur_level - i)/2 + 1] / A[i/2 - 1][0];
        B[i/2][(cur_level - i)/2] =
          A[i/2 - 1][(cur_level - i)/2] / A[i/2 - 1][0];
      }
    }

    // Compute continued fraction coefficients.
    bool reached_zero = false;
    for (uint32_t cur_level = 0; cur_level < 2*num_coeffs; cur_level += 2) {
      if (reached_zero) {
        A[cur_level/2][0] = 0.0;
      }

      if (A[cur_level/2][0] == 0.0) {
        reached_zero = true;
      }

      double pade_coeff;
      pade_coeff = A[cur_level/2][0];

      assert(!std::isnan(pade_coeff) && !std::isinf(pade_coeff));

      coeffs.cont_frac[cur_level/2] = pade_coeff;
    }

    // Compute pade coefficients (as ratio of polynomials).
    std::vector<std::vector<double> > &poly_numerator = coeffs.poly_numerator;
    std::vector<std::vector<double> > &poly_denominator =
      coeffs.poly_denominator;

    poly_numerator[0][0] = A[0][0];
    poly_numerator[1][0] = A[0][0];

    poly_denominator[0][0] = 1.0;
    poly_denominator[1][0] = 1.0;
    poly_denominator[1][1] = A[1][0];

    for (uint32_t k = 2; k < num_coeffs; ++k) {
      if (A[k][0] == 0.0) {
        break;
      }

      assert(k >= 2);
      poly_numerator[k][0] = poly_numerator[k - 1][0];
      poly_denominator[k][0] = poly_denominator[k - 1][0];

      for (uint32_t m = 1; m < num_coeffs; ++m) {
        poly_numerator[k][m] = poly_numerator[k - 1][m]
                             + A[k][0] * poly_numerator[k - 2][m - 1];
        poly_denominator[k][m] = poly_denominator[k - 1][m]
                             + A[k][0] * poly_denominator[k - 2][m - 1];
      }
    }

    return coeffs;
}

PadeRoots ComputeRelevantRootsForConf(uint64_t num_coeffs,
                                      double defect_threshold,
                                      Coeffs const &coeffs) {
  PadeRoots roots;

  // This happens when all the coefficients are 0.
  // In this case, don't try to compute the roots.
  assert(coeffs.cont_frac.size() > 0);
  if (coeffs.cont_frac[0] == 0) {
    return roots;
  }

  bool range_covered = false;
  bool root_is_relevant;

  std::vector<std::vector<double> > const &poly_numerator =
    coeffs.poly_numerator;
  std::vector<std::vector<double> > const &poly_denominator =
    coeffs.poly_denominator;

  std::vector<double> gaps_l, gaps_r, gaps_l_old, gaps_r_old;

  assert(num_coeffs >= 2);
  int k = num_coeffs - 1;

  while (poly_denominator[k][0] == 0.0) {
    --k;
  }
  int first_degree_to_check = k;

  while (!range_covered) {
    assert(k >= 0);

    roots.numerator.push_back(std::vector<double>());
    roots.denominator.push_back(std::vector<double>());

    int degree = (k + 1) / 2;

    gaps_l_old = gaps_l;
    gaps_l.clear();

    gaps_r_old = gaps_r;
    gaps_r.clear();

    double *poly_coeffs_n = new double[(num_coeffs + 1) / 2 + 1];
    double *poly_coeffs_d = new double[(num_coeffs + 1) / 2 + 1];
    double *all_roots_n   = new double[num_coeffs];
    double *all_roots_d   = new double[num_coeffs];

    for (int m = 0; m <= degree; ++m) {
      poly_coeffs_n[m] = poly_numerator[k][degree - m];
    }

    for (int m = 0; m <= degree; ++m) {
      poly_coeffs_d[m] = poly_denominator[k][degree - m];
    }

    gsl_poly_complex_workspace *w =
      gsl_poly_complex_workspace_alloc(degree + 1);
    gsl_poly_complex_solve(poly_coeffs_n, degree + 1, w, all_roots_n);
    gsl_poly_complex_solve(poly_coeffs_d, degree + 1, w, all_roots_d);
    gsl_poly_complex_workspace_free(w);

    for (int m = 0; m < degree; m++) {
      root_is_relevant = false;

      if (all_roots_n[2*m] >= 0 && all_roots_n[2*m + 1] == 0.0) {
        if (k == first_degree_to_check) {
          roots.numerator.back().push_back(all_roots_n[2*m]);
          assert(all_roots_n[2*m] >= 0);

          gaps_l.push_back(std::max(all_roots_n[2*m] - defect_threshold,
                           0.0));
          gaps_r.push_back(all_roots_n[2*m] + defect_threshold);
        } else {
          for (size_t i = 0; i < gaps_l_old.size(); ++i) {
            if (std::abs(all_roots_n[2*m] - gaps_l_old[i])
                  < defect_threshold ||
                std::abs(all_roots_n[2*m] - gaps_r_old[i])
                  < defect_threshold) {
              gaps_l.push_back(
                  std::max(gaps_l_old[i],
                           all_roots_n[2*m] - defect_threshold));
              gaps_r.push_back(
                  std::min(all_roots_n[2*m] + defect_threshold,
                           gaps_r_old[i]));
              if (!root_is_relevant) {
                root_is_relevant = true;
                roots.numerator.back().push_back(all_roots_n[2*m]);
                assert(all_roots_n[2*m] >= 0);
              }
            }
          }
        }
      }
    }

    for (int m = 0; m < degree; ++m) {
      root_is_relevant = false;
      if (all_roots_d[2*m] >= 0 && all_roots_d[2*m + 1] == 0.0) {
        if (k == first_degree_to_check) {
          roots.denominator.back().push_back(all_roots_d[2*m]);
          assert(all_roots_d[2*m] >= 0);

          gaps_l.push_back(
              std::max(all_roots_d[2*m] - defect_threshold, 0.0));
          gaps_r.push_back(all_roots_d[2*m] + defect_threshold);
        } else {
          for (size_t i = 0; i < gaps_l_old.size(); ++i) {
            if (std::abs(all_roots_d[2*m]-gaps_l_old[i]) <
                  defect_threshold ||
                std::abs(all_roots_d[2*m]-gaps_r_old[i]) <
                  defect_threshold) {
              gaps_l.push_back(
                  std::max(gaps_l_old[i],
                           all_roots_d[2*m] - defect_threshold));
              gaps_r.push_back(
                  std::min(all_roots_d[2*m] + defect_threshold,
                           gaps_r_old[i]));
              if (!root_is_relevant) {
                root_is_relevant = true;
                roots.denominator.back().push_back(all_roots_d[2*m]);
                assert(all_roots_d[2*m] >= 0);
              }
            }
          }
        }
      }
    }

    --k;
    range_covered = (gaps_l.size() == 0) || ((k + 1) / 2 == 0);
  }

  // Numerator and denominator should have same size.
  assert(roots.numerator.size() == roots.denominator.size());

  if (roots.numerator.size() > 0 && roots.denominator.size() > 0) {
    // In most cases, the final vectors added to numerator and
    // denominator will have no elements. Remove these vectors for
    // consistency.
    // The only case when the final vectors will be non-empty is when
    // there are relevant roots for all coefficients (except the last
    // (0th) coefficient since that does not have any roots).

    if (roots.numerator.back().size() == 0 &&
        roots.denominator.back().size() == 0) {
      roots.numerator.pop_back();
      roots.denominator.pop_back();
    }
  }

  return roots;
}
