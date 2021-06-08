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

#include "table_gen/linear_solve.h"

#include <stddef.h>
#include <stdint.h>

#include <algorithm>
#include <cmath>

static double const kMaxIterations = 100000;
static double const kTol = 1e-30;

void LinearSolve(std::vector<double> const &Am,
                 std::vector<size_t> const &row_indexes,
                 std::vector<size_t> const &col_indexes,
                 std::vector<double> const &bv,
                 std::vector<double> &x_solution) {
  size_t n = bv.size();
  size_t len_Am = Am.size();

  std::vector<double> x_solution_old(n, 0.0);

  for (uint64_t i = 0; i < kMaxIterations; ++i) {
    uint64_t entry = 0;
    for (uint64_t j = 0; j < n; ++j) {
      double r = 0;
      double s = 0;
      while (entry < len_Am && row_indexes[entry] == j) {
        if (col_indexes[entry] + 1 <= j) {
          r -= Am[entry] * x_solution[col_indexes[entry]];
        } else if (col_indexes[entry] >= j + 1) {
          s -= Am[entry] * x_solution_old[col_indexes[entry]];
        }
        ++entry;
      }
      x_solution[j] = r + s + bv[j];
    }

    double norm = 0;
    for (uint64_t j = 0; j < n; ++j) {
      norm = std::max(norm, std::abs(x_solution[j] - x_solution_old[j]));
    }

    if (norm < kTol) {
      break;
    }

    x_solution.swap(x_solution_old);
  }
}
