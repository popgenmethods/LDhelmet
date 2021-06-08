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

#ifndef LDHELMET_TABLE_GEN_LINEAR_SOLVE_H_
#define LDHELMET_TABLE_GEN_LINEAR_SOLVE_H_

#include <stddef.h>

#include <vector>

// Iterative solver for Ax = b
//   where A is n x n, x is n x 1 and and b is n x 1.

// Linear solver:
// The solution is returned in x.
// Am contains the (non-zero) entries in the matrix.
// row_indexes maps entries to rows.
// col_indexes maps entries to columns.
// n is the number of rows/columns in Am and bv.
// len_Am is the number of entries in Am.
void LinearSolve(std::vector<double> const &Am,
                 std::vector<size_t> const &row_indexes,
                 std::vector<size_t> const &col_indexes,
                 std::vector<double> const &bv,
                 std::vector<double> &x_solution);

#endif  // LDHELMET_TABLE_GEN_LINEAR_SOLVE_H_
