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

#include "pade/compute_g.h"

#include <stdint.h>

#include <vector>

#include "common/conf.h"
#include "pade/one_locus.h"
#include "pade/subtable.h"

double ComputeG(double theta, Table *table,
                uint32_t m, int ind, Conf const &conf) {
  assert(ind < 0 || (static_cast<int>(m) % 2) == (ind % 2));

  uint32_t const kK = 2;
  uint32_t const kL = 2;

  double const theta_a = theta;
  double const theta_b = theta;

  assert(ind >=0);

  assert(static_cast<uint32_t>(conf.ab_) + conf.aB_ + conf.Ab_ + conf.AB_
           == m);

  double ret = 0.0;

  std::vector<double> r_i(kK, 0.0);
  std::vector<double> r_j(kL, 0.0);

  r_i[0] = conf.ab_ + conf.aB_;
  r_i[1] = conf.Ab_ + conf.AB_;

  r_j[0] = conf.ab_ + conf.Ab_;
  r_j[1] = conf.aB_ + conf.AB_;

  uint32_t a = conf.a_ + conf.A_;
  uint32_t b = conf.b_ + conf.B_;

  // Base case: one-locus sampling distribution.
  if (m == 0 && ind == 0) {
    std::vector<int> tmp_a(2);
    tmp_a[0] = conf.a_;
    tmp_a[1] = conf.A_;

    std::vector<int> tmp_b(2);
    tmp_b[0] = conf.b_;
    tmp_b[1] = conf.B_;

    ret = OneLocusSF(theta_a, tmp_a) * OneLocusSF(theta_b, tmp_b);
  } else if (m == 0 && ind > 0) {
    ret = 0.0;
  } else if (m == 1 && (a == 0 || b == 0)) {
    ret = 0.0;
  } else {
    for (uint32_t i = 0; i < kK; ++i) {
      for (uint32_t j = 0; j < kL; ++j) {
        double r_ij = static_cast<double>(conf.Getr(i, j));

        if (m >= 2 && conf.Getr(i, j) >= 2) {
          ret += r_ij * (r_ij - 1.0)
               * GetTable(*table, m-2, ind,
                          conf.Adda(i, 1)
                              .Addb(j, 1)
                              .Addr(i, j, -2));
        }
        for (uint32_t l = 0; l < kL; ++l) {
          double r_il = static_cast<double>(conf.Getr(i, l));

          if (m >= 2 &&
             (j == l ? conf.Getr(i, j) >= 2
                     : conf.Getr(i, j) >= 1 && conf.Getr(i, l) >= 1)) {
            ret -= r_ij * (r_il - delta(j, l))
                 * GetTable(*table, m - 2, ind,
                            conf.Adda(i, 1)
                                .Addb(j, 1)
                                .Addb(l, 1)
                                .Addr(i, j, -1)
                                .Addr(i, l, -1));
          }
        }

        for (uint32_t k = 0; k < kK; ++k) {
          double r_kj = static_cast<double>(conf.Getr(k, j));

          if (m >= 2 &&
              (i == k ? conf.Getr(i, j) >= 2
                      : conf.Getr(i, j) >= 1 && conf.Getr(k, j) >= 1)) {
            ret -= r_ij * (r_kj - delta(i, k))
                 * GetTable(*table, m - 2, ind,
                            conf.Adda(i, 1)
                                .Adda(k, 1)
                                .Addb(j, 1)
                                .Addr(i, j, -1)
                                .Addr(k, j, -1));
          }
        }

        for (uint32_t k = 0; k < kK; ++k) {
          for (uint32_t l = 0; l < kL; ++l) {
            double r_kl = static_cast<double>(conf.Getr(k, l));

            if (m >= 2 &&
                (i == k && j == l ?
                   conf.Getr(i, j) >= 2
                 : conf.Getr(i, j) >= 1 && conf.Getr(k, l) >= 1)) {
              ret += r_ij * (r_kl - delta(i, k) * delta(j, l))
                   * GetTable(*table, m - 2, ind,
                              conf.Adda(i, 1)
                                  .Adda(k, 1)
                                  .Addb(j, 1)
                                  .Addb(l, 1)
                                  .Addr(i, j, -1)
                                  .Addr(k, l, -1));
            }
          }
        }
      }
    }

    for (uint32_t i = 0; i < kK; ++i) {
        for (uint32_t j = 0; j < kL; ++j) {
            double r_ij = static_cast<double>(conf.Getr(i, j));

            if (m >= 1 && ind >= 1 && conf.Getr(i, j) >= 1) {
              ret += r_ij * (r_ij - 1.0)
                 * GetTable(*table, m - 1, ind - 1,
                            conf.Addr(i, j, -1));
            }

            if (m >= 1 && ind >= 1 && conf.Getr(i, j) >= 1) {
              ret -= 2.0 * r_ij * (r_i[i] - 1.0)
                   * GetTable(*table, m - 1, ind - 1,
                              conf.Addb(j, 1)
                                  .Addr(i, j, -1));
            }

            if (m >= 1 && ind >= 1 && conf.Getr(i, j) >= 1) {
              ret -= 2.0 * r_ij * (r_j[j] - 1.0)
                   * GetTable(*table, m - 1, ind - 1,
                              conf.Adda(i, 1)
                                  .Addr(i, j, -1));
            }

            for (uint32_t k = 0; k < kK; ++k) {
              for (uint32_t l = 0; l < kL; ++l) {
                double r_kj = static_cast<double>(conf.Getr(k, j));
                double r_il = static_cast<double>(conf.Getr(i, l));

                bool valid;
                if (k == i) {
                  if (j == l) {
                    valid = conf.Getr(i, j) >= 1;
                  } else {
                    valid = conf.Getr(i, l) >= 1;
                  }
                } else {
                  if (j == l) {
                    valid = conf.Getr(k, j) >= 1;
                  } else {
                    valid = conf.Getr(k, j) >= 1 && conf.Getr(i, l) >= 1;
                  }
                }

                if (m >= 1 && ind >= 1 && valid) {
                  ret += 2.0 * r_kj * (r_il - delta(i, k) * delta(j, l))
                       * GetTable(*table, m-1, ind-1,
                                  conf.Adda(k, 1)
                                      .Addb(l, 1)
                                      .Addr(k, j, -1)
                                      .Addr(i, l, -1)
                                      .Addr(i, j, 1));
                }
              }
            }

            if (m >= 1 && ind >= 1 && conf.Getr(i, j) >= 1) {
              ret += 2.0 * (m - 1) * r_ij
                   * GetTable(*table, m - 1, ind - 1,
                              conf.Adda(i, 1)
                                  .Addb(j, 1)
                                  .Addr(i, j, -1));
            }
        }
    }

    for (uint32_t i = 0; i < kK; ++i) {
      double a_i = static_cast<double>(conf.Geta(i));
      if (ind >= 2 && conf.Geta(i) >= 1) {
        ret += a_i * (a_i + 2.0 * r_i[i] - 1.0)
             * GetTable(*table, m, ind - 2,
                        conf.Adda(i, -1));
      }
    }

    for (uint32_t j = 0; j < kL; ++j) {
      double b_j = static_cast<double>(conf.Getb(j));
      if (ind >= 2 && conf.Getb(j) >= 1) {
        ret += b_j * (b_j + 2.0 * r_j[j] - 1.0)
             * GetTable(*table, m, ind - 2,
                        conf.Addb(j, -1));
      }
    }

    for (uint32_t i = 0; i < kK; ++i) {
      double a_i = static_cast<double>(conf.Geta(i));
      for (uint32_t j = 0; j < kL; ++j) {
        double b_j = static_cast<double>(conf.Getb(j));
        for (uint32_t k = 0; k < kK; ++k) {
          double r_kj = static_cast<double>(conf.Getr(k, j));

          if (ind >= 2 && conf.Geta(i) >= 1 &&
              (i == k ? true : conf.Getr(k, j) >= 1)) {
            ret -= 2.0 * a_i * r_kj
                 * GetTable(*table, m, ind - 2,
                            conf.Adda(i, -1)
                                .Adda(k, 1)
                                .Addr(k, j, -1)
                                .Addr(i, j, 1));
          }
        }
        for (uint32_t l = 0; l < kL; ++l) {
          double r_il = static_cast<double>(conf.Getr(i, l));

          if (ind >= 2 && conf.Getb(j) >= 1 &&
              (l == j ? true : conf.Getr(i, l) >= 1)) {
            ret -= 2.0 * b_j * r_il
                 * GetTable(*table, m, ind - 2,
                            conf.Addb(j, -1)
                                .Addb(l, 1)
                                .Addr(i, l, -1)
                                .Addr(i, j, 1));
          }
        }
      }
    }

    if (ind >= 2 && conf.Geta(1) == 1 && r_i[1] == 0) {
      ret += theta_a * GetTable(*table, m, ind - 2,
                                conf.Adda(1, -1)
                                    .Adda(0, 1));
    } else if (conf.Geta(1) == 0 && r_i[1] == 1) {
      for (uint32_t j = 0; j < kL; ++j) {
        if (ind >= 2 && conf.Getr(1, j) >= 1) {
          ret += theta_a * conf.Getr(1, j)
               * GetTable(*table, m, ind - 2,
                          conf.Addr(1, j, -1)
                              .Addr(0, j, 1));
        }
      }
    }

    if (ind >= 2 && conf.Getb(1) == 1 && r_j[1] == 0) {
      ret += theta_b * GetTable(*table, m, ind - 2,
                                conf.Addb(1, -1)
                                    .Addb(0, 1));
    } else if (conf.Getb(1) == 0 && r_j[1] == 1) {
      for (uint32_t i = 0; i < kK; ++i) {
        if (ind >= 2 && conf.Getr(i, 1) >= 1) {
          ret += theta_b * conf.Getr(i, 1)
               * GetTable(*table, m, ind - 2,
                          conf.Addr(i, 1, -1)
                              .Addr(i, 0, 1));
        }
      }
    }

    if (ind >= 2) {
      ret -= ((a + m) * (a + m + theta_a - 1.0)
               + (b + m) * (b + m + theta_b - 1.0)
               - m * (m - 3.0))
           * GetTable(*table, m, ind - 2, conf);
    }

    for (uint32_t i = 0; i < kK; ++i) {
      double a_i = static_cast<double>(conf.Geta(i));
      for (uint32_t j = 0; j < kL; ++j) {
        double b_j = static_cast<double>(conf.Getb(j));
        if (ind >= 3 && conf.Geta(i) >= 1 && conf.Getb(j) >= 1) {
          ret += 2.0 * a_i * b_j
               * GetTable(*table, m + 1, ind - 3,
                          conf.Adda(i, -1)
                              .Addb(j, -1)
                              .Addr(i, j, 1));
        }
      }
    }

    ret /= m;
  }

  SetTable(table, m, ind, conf, ret);

  return ret;
}
