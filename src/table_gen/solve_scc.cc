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

#include "table_gen/solve_scc.h"

#include <algorithm>
#include <vector>

#include "common/conf.h"
#include "common/ncr.h"
#include "common/vector_definitions.h"
#include "table_gen/linear_solve.h"
#include "table_gen/table_management.h"

uint64_t ComputeNumConfs(uint32_t degree) {
  uint64_t num = 0;
  for (uint64_t w = 0; w < degree / 2 + 1; ++w) {
    num += Ncr(static_cast<uint64_t>(degree) - 2 * w +4 - 1, 4 - 1)
             * Ncr(w + 4 - 1, 4 - 1);
  }
  return num;
}

uint64_t ComputeNumConfsSCC(uint32_t a_mar,
                            uint32_t A_mar,
                            uint32_t b_mar,
                            uint32_t B_mar) {
  uint64_t num_confs_scc = 0;
  uint32_t min_ab = std::min(a_mar, b_mar);
  for (uint32_t ab = 0; ab <= min_ab; ++ab) {
    uint32_t min_aB = std::min(a_mar-ab, B_mar);
    for (uint32_t aB = 0; aB <= min_aB; ++aB) {
      uint32_t min_Ab = std::min(b_mar - ab, A_mar);
      for (uint32_t Ab = 0; Ab <= min_Ab; ++Ab) {
        uint32_t min_AB = std::min(A_mar - Ab, B_mar - aB);
        num_confs_scc += min_AB + 1;
      }
    }
  }
  return num_confs_scc;
}

IndexTables ComputeIndexTables(uint32_t a_mar,
                               uint32_t A_mar,
                               uint32_t b_mar,
                               uint32_t B_mar) {
  IndexTables index_tables;

  // Compute table_Ab_.
  index_tables.table_Ab_.resize(std::min(a_mar, B_mar) + 1);
  for (uint32_t aB = 0; aB < index_tables.table_Ab_.size(); ++aB) {
    index_tables.table_Ab_[aB].resize(std::min(b_mar, A_mar) + 2);
    size_t sum = 0;
    for (uint32_t Ab = 0;
         Ab < index_tables.table_Ab_[aB].size() - 1;
         ++Ab) {
      index_tables.table_Ab_[aB][Ab] = sum;
      sum += std::min(A_mar - Ab, B_mar - aB) + 1;
    }
    index_tables.table_Ab_[aB][index_tables.table_Ab_[aB].size() - 1] = sum;
  }

  // Compute table_aB_.
  index_tables.table_aB_.resize(std::min(a_mar, b_mar) + 1);
  for (uint32_t ab = 0; ab < index_tables.table_aB_.size(); ++ab) {
    index_tables.table_aB_[ab].resize(std::min(a_mar, B_mar) + 2);
    size_t sum = 0;
    for (uint32_t aB = 0; aB < index_tables.table_aB_[ab].size() - 1; ++aB) {
        index_tables.table_aB_[ab][aB] = sum;
        sum += index_tables.table_Ab_[aB][std::min(b_mar - ab, A_mar) + 1];
    }
    index_tables.table_aB_[ab][index_tables.table_aB_[ab].size() - 1] = sum;
  }

  // Compute table_ab_.
  index_tables.table_ab_.resize(std::min(a_mar, b_mar) + 2);
  {
    size_t sum = 0;
    for (uint32_t ab = 0; ab < index_tables.table_ab_.size() - 1; ++ab) {
      index_tables.table_ab_[ab] = sum;
      sum += index_tables.table_aB_[ab][std::min(a_mar - ab, B_mar) + 1];
    }
    index_tables.table_ab_[index_tables.table_ab_.size() - 1] = sum;
  }

  return index_tables;
}

// Configuration to index within SCC.
size_t ConfToIndex(IndexTables const &index_tables, Conf const &conf) {
  assert(conf.ab_ < index_tables.table_ab_.size());
  assert(conf.ab_ < index_tables.table_aB_.size());
  assert(conf.aB_ < index_tables.table_aB_[conf.ab_].size());
  assert(conf.aB_ < index_tables.table_Ab_.size());
  assert(conf.Ab_ < index_tables.table_Ab_[conf.aB_].size());

  return index_tables.table_ab_[conf.ab_]
       + index_tables.table_aB_[conf.ab_][conf.aB_]
       + index_tables.table_Ab_[conf.aB_][conf.Ab_]
       + conf.AB_;
}

Conf IndexToConf(IndexTables const &index_tables, size_t index,
                 uint32_t a_mar, uint32_t A_mar,
                 uint32_t b_mar, uint32_t B_mar) {
  Conf conf;

  // Set ab.
  while (index >= index_tables.table_ab_[conf.ab_ + 1]) {
    ++conf.ab_;
  }
  index -= index_tables.table_ab_[conf.ab_];

  assert(conf.ab_ != index_tables.table_ab_.size() - 1);
  assert(conf.ab_ <= std::min(a_mar, b_mar));

  // Set aB.
  while (index >= index_tables.table_aB_[conf.ab_][conf.aB_ + 1]) {
    ++conf.aB_;
  }
  index -= index_tables.table_aB_[conf.ab_][conf.aB_];

  assert(conf.aB_ != index_tables.table_aB_[conf.ab_].size() - 1);
  assert(conf.aB_ <= std::min(a_mar-conf.ab_, B_mar));

  // Set Ab.
  while (index >= index_tables.table_Ab_[conf.aB_][conf.Ab_ + 1]) {
    ++conf.Ab_;
  }
  index -= index_tables.table_Ab_[conf.aB_][conf.Ab_];

  assert(conf.Ab_ != index_tables.table_Ab_[conf.aB_].size() - 1);
  assert(conf.Ab_ <= std::min(b_mar - conf.ab_, A_mar));

  // Set AB.
  conf.AB_ = index;
  assert(conf.AB_ <= std::min(A_mar - conf.Ab_, B_mar - conf.aB_));
  // Set a.
  conf.a_ = a_mar - conf.ab_ - conf.aB_;
  // Set A.
  conf.A_ = A_mar - conf.Ab_ - conf.AB_;
  // Set b.
  conf.b_ = b_mar - conf.ab_ - conf.Ab_;
  // Set B.
  conf.B_ = B_mar - conf.aB_ - conf.AB_;

  return conf;
}

void AdjustAm(IndexTables const &index_tables,
              size_t conf_index,
              double denominator,
              double coeff,
              Conf const &tmp_conf,
              std::vector<double> &Am,
              std::vector<size_t> &row_indexes,
              std::vector<size_t> &col_indexes) {
  // Conf to column index.
  size_t tmp_conf_index = ConfToIndex(index_tables, tmp_conf);

  Am.push_back(-1.0 * (coeff / denominator));
  row_indexes[Am.size() - 1] = conf_index;
  col_indexes[Am.size() - 1] = tmp_conf_index;
}

void ConstructSystem(Vec8 const &table,
                     IndexTables const &index_tables,
                     double theta,
                     double rho,
                     uint32_t a_mar,
                     uint32_t A_mar,
                     uint32_t b_mar,
                     uint32_t B_mar,
                     std::vector<double> &Am,
                     std::vector<size_t> &row_indexes,
                     std::vector<size_t> &col_indexes,
                     std::vector<double> &bv) {
  Conf conf;
  uint32_t min_ab = std::min(a_mar, b_mar);
  for (conf.ab_ = 0; conf.ab_ <= min_ab; ++conf.ab_) {
    uint32_t min_aB = std::min(a_mar - conf.ab_, B_mar);
    for (conf.aB_ = 0; conf.aB_ <= min_aB; ++conf.aB_) {
      uint32_t min_Ab = std::min(b_mar - conf.ab_, A_mar);
      for (conf.Ab_= 0; conf.Ab_ <= min_Ab; ++conf.Ab_) {
        uint32_t min_AB = std::min(A_mar - conf.Ab_, B_mar - conf.aB_);
        for (conf.AB_ = 0; conf.AB_ <= min_AB; ++conf.AB_) {
          conf.a_ = a_mar - conf.ab_ - conf.aB_;
          conf.A_ = A_mar - conf.Ab_ - conf.AB_;
          conf.b_ = b_mar - conf.ab_ - conf.Ab_;
          conf.B_ = B_mar - conf.aB_ - conf.AB_;

          uint32_t degree = a_mar + A_mar + b_mar + B_mar;
          assert((static_cast<uint32_t>(conf.a_) +
                   conf.ab_ + conf.aB_) == a_mar);

          assert((static_cast<uint32_t>(conf.A_) +
                   conf.Ab_ + conf.AB_) == A_mar);

          assert((static_cast<uint32_t>(conf.b_) +
                   conf.ab_ + conf.Ab_) == b_mar);

          assert((static_cast<uint32_t>(conf.B_) +
                   conf.aB_ + conf.AB_) == B_mar);

          assert((static_cast<uint32_t>(conf.a_) +
                  conf.A_+conf.b_+conf.B_)
                 + 2*(conf.ab_ + conf.aB_ + conf.Ab_ + conf.AB_)
                   == degree);

          assert(a_mar + A_mar + b_mar + B_mar == degree);

          // conf_index is assumed to increment consecutively with
          // each iteration.
          size_t conf_index = ConfToIndex(index_tables, conf);

          // Fill out b vector and A matrix.
          uint32_t a_tot = conf.a_ + conf.A_;  //  a + A
          uint32_t b_tot = conf.b_ + conf.B_;  //  b + B
          uint32_t c_a = conf.ab_ + conf.aB_;  //  ab + aB
          uint32_t c_A = conf.Ab_ + conf.AB_;  //  Ab + AB
          uint32_t c_b = conf.ab_ + conf.Ab_;  //  ab + Ab
          uint32_t c_B = conf.aB_ + conf.AB_;  //  aB + AB

          // ab + aB + Ab + AB
          uint32_t c_tot = conf.ab_ + conf.aB_ + conf.Ab_ + conf.AB_;

          // a + A + b + B + ab + aB + Ab + AB
          uint32_t n = a_tot + b_tot + c_tot;

          double denominator = static_cast<double>(n) * static_cast<double>(n-1)
                             + theta * static_cast<double>(a_tot + c_tot)
                             + theta * static_cast<double>(b_tot + c_tot)
                             + rho * static_cast<double>(c_tot);

          // Modify bv.
          bv[conf_index] = 0.0;

          // *** Transitions that depend on degree-1 ***
          // Transitions I (coalescence of haps defined only at locus 1).
          if (conf.a_ > 0 && conf.a_ - 1 + 2*c_a > 0) {
            AdjustBv(table,
                     conf_index,
                     denominator,
                     static_cast<double>(conf.a_ * (conf.a_ - 1 + 2*c_a)),
                     Conf(conf.a_ - 1,
                          conf.A_,
                          conf.b_,
                          conf.B_,
                          conf.ab_,
                          conf.aB_,
                          conf.Ab_,
                          conf.AB_),
                     bv);
          }

          if (conf.A_ > 0 && conf.A_ - 1 + 2*c_A > 0) {
            AdjustBv(table,
                     conf_index,
                     denominator,
                     static_cast<double>(conf.A_ * (conf.A_ - 1 + 2*c_A)),
                     Conf(conf.a_,
                          conf.A_ - 1,
                          conf.b_,
                          conf.B_,
                          conf.ab_,
                          conf.aB_,
                          conf.Ab_,
                          conf.AB_),
                     bv);
          }

          // Transitions II (coalescence of haps defined only at locus 2).
          if (conf.b_ > 0 && conf.b_ - 1 + 2*c_b > 0) {
              AdjustBv(table,
                       conf_index,
                       denominator,
                       static_cast<double>(conf.b_ * (conf.b_ - 1 + 2*c_b)),
                       Conf(conf.a_,
                            conf.A_,
                            conf.b_ - 1,
                            conf.B_,
                            conf.ab_,
                            conf.aB_,
                            conf.Ab_,
                            conf.AB_),
                       bv);
          }
          if (conf.B_ > 0 && conf.B_ - 1 + 2*c_B > 0) {
            AdjustBv(table,
                     conf_index,
                     denominator,
                     static_cast<double>(conf.B_ * (conf.B_ - 1 + 2*c_B)),
                     Conf(conf.a_,
                          conf.A_,
                          conf.b_,
                          conf.B_ - 1,
                          conf.ab_,
                          conf.aB_,
                          conf.Ab_,
                          conf.AB_),
                     bv);
          }

          // Transitions III (coalescence of fully defined haps).
          if (conf.ab_ > 1) {
            AdjustBv(table,
                     conf_index,
                     denominator,
                     static_cast<double>(conf.ab_ * (conf.ab_ - 1)),
                     Conf(conf.a_,
                          conf.A_,
                          conf.b_,
                          conf.B_,
                          conf.ab_ - 1,
                          conf.aB_,
                          conf.Ab_,
                          conf.AB_),
                     bv);
          }

          if (conf.aB_ > 1) {
            AdjustBv(table,
                     conf_index,
                     denominator,
                     static_cast<double>(conf.aB_  *(conf.aB_ - 1)),
                     Conf(conf.a_,
                          conf.A_,
                          conf.b_,
                          conf.B_,
                          conf.ab_,
                          conf.aB_ - 1,
                          conf.Ab_,
                          conf.AB_),
                     bv);
          }

          if (conf.Ab_ > 1) {
            AdjustBv(table,
                     conf_index,
                     denominator,
                     static_cast<double>(conf.Ab_*(conf.Ab_ - 1)),
                     Conf(conf.a_,
                          conf.A_,
                          conf.b_,
                          conf.B_,
                          conf.ab_,
                          conf.aB_,
                          conf.Ab_ - 1,
                          conf.AB_),
                     bv);
          }

          if (conf.AB_ > 1) {
            AdjustBv(table,
                     conf_index,
                     denominator,
                     static_cast<double>(conf.AB_*(conf.AB_ - 1)),
                     Conf(conf.a_,
                          conf.A_,
                          conf.b_,
                          conf.B_,
                          conf.ab_,
                          conf.aB_,
                          conf.Ab_,
                          conf.AB_ - 1),
                     bv);
          }

          // *** Transitions that depend on current degree in SCCs
          //     that should be already solved ***
          // Transitions V, VI, VII, VIII (mutation).
          // Mutation only occurs when (A+Ab+AB==1 || B+aB+AB==1),
          // in which case, the mutation changes the derived allele
          // to the ancestral allele.
          // A configuration may allow no mutations, 1 mutation, or 2
          // mutations.
          // Note that the probabilities for these mutations should
          // be known already because they should have been solved
          // for previously.
          if (A_mar == 1) {  // Mutation at locus 1.
            if (conf.A_ == 1) {
              AdjustBv(table,
                       conf_index,
                       denominator,
                       theta * static_cast<double>(conf.A_),
                       Conf(conf.a_ + 1,
                            conf.A_ - 1,
                            conf.b_,
                            conf.B_,
                            conf.ab_,
                            conf.aB_,
                            conf.Ab_,
                            conf.AB_),
                       bv);
            }

            if (conf.Ab_ == 1) {
              AdjustBv(table,
                       conf_index,
                       denominator,
                       theta * static_cast<double>(conf.Ab_),
                       Conf(conf.a_,
                            conf.A_,
                            conf.b_,
                            conf.B_,
                            conf.ab_ + 1,
                            conf.aB_,
                            conf.Ab_ - 1,
                            conf.AB_),
                       bv);
            }

            if (conf.AB_ == 1) {
              AdjustBv(table,
                       conf_index,
                       denominator,
                       theta * static_cast<double>(conf.AB_),
                       Conf(conf.a_,
                            conf.A_,
                            conf.b_,
                            conf.B_,
                            conf.ab_,
                            conf.aB_ + 1,
                            conf.Ab_,
                            conf.AB_ - 1),
                       bv);
            }
          }

          if (B_mar == 1) {  // Mutation at locus 2.
            if (conf.B_ == 1) {
              AdjustBv(table,
                       conf_index,
                       denominator,
                       theta * static_cast<double>(conf.B_),
                       Conf(conf.a_,
                            conf.A_,
                            conf.b_ + 1,
                            conf.B_ - 1,
                            conf.ab_,
                            conf.aB_,
                            conf.Ab_,
                            conf.AB_),
                       bv);
            }
            if (conf.aB_ == 1) {
              AdjustBv(table,
                       conf_index,
                       denominator,
                       theta * static_cast<double>(conf.aB_),
                       Conf(conf.a_,
                            conf.A_,
                            conf.b_,
                            conf.B_,
                            conf.ab_ + 1,
                            conf.aB_ - 1,
                            conf.Ab_,
                            conf.AB_),
                       bv);
            }
            if (conf.AB_ == 1) {
              AdjustBv(table,
                       conf_index,
                       denominator,
                       theta * static_cast<double>(conf.AB_),
                       Conf(conf.a_,
                            conf.A_,
                            conf.b_,
                            conf.B_,
                            conf.ab_,
                            conf.aB_,
                            conf.Ab_ + 1,
                            conf.AB_ - 1),
                       bv);
            }
          }

          // *** Transitions that depend on current degree ***
          // Transitions IV (coalescence of 2 half-defined haps into
          // a fully defined hap).
          if (conf.a_ > 0 && conf.b_ > 0) {
            AdjustAm(index_tables,
                     conf_index,
                     denominator,
                     static_cast<double>(2 * conf.a_ * conf.b_),
                     Conf(conf.a_ - 1,
                          conf.A_,
                          conf.b_ - 1,
                          conf.B_,
                          conf.ab_ + 1,
                          conf.aB_,
                          conf.Ab_,
                          conf.AB_),
                     Am,
                     row_indexes,
                     col_indexes);
          }

          if (conf.a_ > 0 && conf.B_ > 0) {
            AdjustAm(index_tables,
                     conf_index,
                     denominator,
                     static_cast<double>(2 * conf.a_ * conf.B_),
                     Conf(conf.a_ - 1,
                          conf.A_,
                          conf.b_,
                          conf.B_ - 1,
                          conf.ab_,
                          conf.aB_ + 1,
                          conf.Ab_,
                          conf.AB_),
                     Am,
                     row_indexes,
                     col_indexes);
          }

          if (conf.A_ > 0 && conf.b_ > 0) {
            AdjustAm(index_tables,
                     conf_index,
                     denominator,
                     static_cast<double>(2 * conf.A_ * conf.b_),
                     Conf(conf.a_,
                          conf.A_ - 1,
                          conf.b_ - 1,
                          conf.B_,
                          conf.ab_,
                          conf.aB_,
                          conf.Ab_ + 1,
                          conf.AB_),
                     Am,
                     row_indexes,
                     col_indexes);
          }

          if (conf.A_ > 0 && conf.B_ > 0) {
            AdjustAm(index_tables,
                     conf_index,
                     denominator,
                     static_cast<double>(2 * conf.A_ * conf.B_),
                     Conf(conf.a_,
                          conf.A_ - 1,
                          conf.b_,
                          conf.B_ - 1,
                          conf.ab_,
                          conf.aB_,
                          conf.Ab_,
                          conf.AB_ + 1),
                     Am,
                     row_indexes,
                     col_indexes);
          }

          // Transitions IX (recombination).
          if (conf.ab_ > 0) {
            AdjustAm(index_tables,
                     conf_index,
                     denominator,
                     rho * static_cast<double>(conf.ab_),
                     Conf(conf.a_ + 1,
                          conf.A_,
                          conf.b_ + 1,
                          conf.B_,
                          conf.ab_ - 1,
                          conf.aB_,
                          conf.Ab_,
                          conf.AB_),
                     Am,
                     row_indexes,
                     col_indexes);
          }

          if (conf.aB_ > 0) {
            AdjustAm(index_tables,
                     conf_index,
                     denominator,
                     rho * static_cast<double>(conf.aB_),
                     Conf(conf.a_ + 1,
                          conf.A_,
                          conf.b_,
                          conf.B_ + 1,
                          conf.ab_,
                          conf.aB_ - 1,
                          conf.Ab_,
                          conf.AB_),
                     Am,
                     row_indexes,
                     col_indexes);
          }

          if (conf.Ab_ > 0) {
            AdjustAm(index_tables,
                     conf_index,
                     denominator,
                     rho * static_cast<double>(conf.Ab_),
                     Conf(conf.a_,
                          conf.A_ + 1,
                          conf.b_ + 1,
                          conf.B_,
                          conf.ab_,
                          conf.aB_,
                          conf.Ab_ - 1,
                          conf.AB_),
                     Am,
                     row_indexes,
                     col_indexes);
          }
          if (conf.AB_ > 0) {
            AdjustAm(index_tables,
                     conf_index,
                     denominator,
                     rho * static_cast<double>(conf.AB_),
                     Conf(conf.a_,
                          conf.A_ + 1,
                          conf.b_,
                          conf.B_ + 1,
                          conf.ab_,
                          conf.aB_,
                          conf.Ab_,
                          conf.AB_ - 1),
                     Am,
                     row_indexes,
                     col_indexes);
          }

          Am.push_back(1.0);
          col_indexes[Am.size() - 1] = conf_index;
          row_indexes[Am.size() - 1] = conf_index;
        }
      }
    }
  }
}

void SolveSCC(Vec8 *table,
              double theta,
              double rho,
              uint32_t a_mar,
              uint32_t A_mar,
              uint32_t b_mar,
              uint32_t B_mar) {
  // Compute index tables.
  IndexTables index_tables = ComputeIndexTables(a_mar, A_mar, b_mar, B_mar);

  // A matrix:
  //   Each row corresponds to a configuration in the SCC.
  //   Each column corresponds to the coefficient for another configuration
  //   in the SCC according to the recursion equation.
  // b vector:
  //   Each component corresponds to a configuration in the SCC.

  uint64_t num_confs_scc = ComputeNumConfsSCC(a_mar, A_mar, b_mar, B_mar);

  std::vector<double> Am;
  Am.reserve(9 * num_confs_scc);

  std::vector<size_t> row_indexes(9 * num_confs_scc, 0.0);
  std::vector<size_t> col_indexes(9 * num_confs_scc, 0.0);
  std::vector<double> bv(num_confs_scc);

  ConstructSystem(*table,
                  index_tables,
                  theta,
                  rho,
                  a_mar,
                  A_mar,
                  b_mar,
                  B_mar,
                  Am,
                  row_indexes,
                  col_indexes,
                  bv);

  // *** Solve system ***
  std::vector<double> x_solution(num_confs_scc);
  LinearSolve(Am, row_indexes, col_indexes, bv, x_solution);

  // Add results to table.
  for (size_t conf_index = 0; conf_index < num_confs_scc; ++conf_index) {
    Conf conf = IndexToConf(index_tables, conf_index,
                            a_mar, A_mar, b_mar, B_mar);
    SetTable(table, conf, x_solution[conf_index]);
  }
}
