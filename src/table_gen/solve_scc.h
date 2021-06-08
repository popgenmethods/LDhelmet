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

#ifndef LDHELMET_TABLE_GEN_SOLVE_SCC_H_
#define LDHELMET_TABLE_GEN_SOLVE_SCC_H_

#include <stddef.h>
#include <stdint.h>

#include <vector>

#include "common/conf.h"
#include "common/ncr.h"
#include "common/vector_definitions.h"
#include "table_gen/table_management.h"

// IndexTables is used to map from conf to index and from index to conf.
class IndexTables {
 public:
  std::vector<std::vector<size_t> > table_Ab_;
  std::vector<std::vector<size_t> > table_aB_;
  std::vector<size_t> table_ab_;
};

// Compute number of configurations in a given degree.
uint64_t ComputeNumConfs(uint32_t degree);

// Compute number of configurations in a given SCC.
uint64_t ComputeNumConfsSCC(uint32_t a_mar,
                            uint32_t A_mar,
                            uint32_t b_mar,
                            uint32_t B_mar);

// Computes index tables for mapping from conf to index and from index to conf.
IndexTables ComputeIndexTables(uint32_t a_mar,
                               uint32_t A_mar,
                               uint32_t b_mar,
                               uint32_t B_mar);

// Configuration to index within SCC.
size_t ConfToIndex(IndexTables const &index_tables, Conf const &conf);

// Index to configuration within SCC.
Conf IndexToConf(IndexTables const &index_tables, size_t index,
                 uint32_t a_mar, uint32_t A_mar,
                 uint32_t b_mar, uint32_t B_mar);

inline void AdjustBv(Vec8 const &table,
                   size_t conf_index,
                   double denominator,
                   double coeff,
                   Conf const &tmp_conf,
                   std::vector<double> &bv) {
  bv[conf_index] += (coeff / denominator) * GetTable(table, tmp_conf);
}

void AdjustAm(IndexTables const &index_tables,
              size_t conf_index,
              double denominator,
              double coeff,
              Conf const &tmp_conf,
              std::vector<double> &Am,
              std::vector<size_t> &row_indexes,
              std::vector<size_t> &col_indexes);

// Set up system of linear equations for a given configuration.
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
                     std::vector<double> &bv);

// Solve SCC for inductive case.
// Assumes that all SCCs it depends on have already been computed.
void SolveSCC(Vec8 *table,
              double theta,
              double rho,
              uint32_t a_mar,
              uint32_t A_mar,
              uint32_t b_mar,
              uint32_t B_mar);

#endif  // LDHELMET_TABLE_GEN_SOLVE_SCC_H_
