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

#include "table_gen/table_management.h"

#include <stdint.h>

#include <algorithm>

#include "common/conf.h"
#include "common/vector_definitions.h"

uint64_t AllocateMemoryToTable(Vec8 *table, uint32_t degree) {
  // Allocate memory to table for this degree.
  // table[degree][a_mar][A_mar][b_mar][ab][aB][Ab][AB]
  uint64_t num_confs = 0;
  // Size is degree + 1 for index a_mar.
  (*table)[degree] = Vec7(degree + 1, Vec6());
  for (uint32_t a_mar = 0; a_mar <= degree; ++a_mar) {
    // Size is degree - a_mar + 1 for index A_mar.
    (*table)[degree][a_mar].resize(degree - a_mar + 1);
    uint32_t A_mar_max = degree - a_mar;
    for (uint32_t A_mar = 0; A_mar <= A_mar_max; ++A_mar) {
      // Size is degree-a_mar-A_mar+1 for index b_mar.
      (*table)[degree][a_mar][A_mar].resize(degree - a_mar - A_mar + 1);
      // Restrict b_mar to be less than or equal to a_mar to maintain
      // canonical form.
      uint32_t b_mar_max = std::min(degree - a_mar - A_mar, a_mar);
      for (uint32_t b_mar = 0; b_mar <= b_mar_max; ++b_mar) {
        uint32_t B_mar = degree - a_mar - A_mar - b_mar;
        // Skip this SCC because we want (a_mar, A_mar) >= (b_mar, B_mar)
        // (in lexicographic order).
        if (b_mar == a_mar && B_mar > A_mar) {
          continue;
        }

        // At this point, a_mar, A_mar, b_mar, and B_mar are fixed.
        // Enumerate all configurations consistent with this SCC.
        uint32_t min_ab = std::min(a_mar, b_mar);

        // Size min(a_mar,b_mar)+1 for index ab.
        (*table)[degree][a_mar][A_mar][b_mar].resize(min_ab + 1);
        for (uint32_t ab = 0; ab <= min_ab; ++ab) {
          uint32_t min_aB = std::min(a_mar - ab, B_mar);
          // Size is min(a_mar - ab, B_mar) + 1 for index aB.
          (*table)[degree][a_mar][A_mar][b_mar][ab].resize(min_aB + 1);
          for (uint32_t aB = 0; aB <= min_aB; ++aB) {
            uint32_t min_Ab = std::min(b_mar - ab, A_mar);
            // Size is min(min(b_mar - ab, A_mar)+1, aB) for index Ab
            // (with symmetry optimization).
            (*table)[degree][a_mar][A_mar][b_mar][ab][aB].resize(min_Ab + 1);
            for (uint32_t Ab = 0; Ab <= min_Ab; ++Ab) {
              uint32_t min_AB = std::min(A_mar - Ab, B_mar - aB);
              // Size is min(A_mar-Ab,B_mar-aB)+1 for index AB.
              (*table)[degree][a_mar][A_mar]
                              [b_mar][ab][aB][Ab].resize(min_AB + 1);
              num_confs += min_AB + 1;
            }
          }
        }
      }
    }
  }
  return num_confs;
}

Conf GetCanon(Conf conf) {
  uint32_t a_mar = conf.a_ + conf.ab_ + conf.aB_;
  uint32_t A_mar = conf.A_ + conf.Ab_ + conf.AB_;
  uint32_t b_mar = conf.b_ + conf.ab_ + conf.Ab_;
  uint32_t B_mar = conf.B_ + conf.aB_ + conf.AB_;
  if (!(a_mar == b_mar ? A_mar >= B_mar : a_mar > b_mar)) {
    std::swap(conf.aB_, conf.Ab_);
    std::swap(conf.a_, conf.b_);
    std::swap(conf.A_, conf.B_);
  }
  return conf;
}

double GetTable(Vec8 const &table, Conf const &in_conf) {
  Conf conf = GetCanon(in_conf);
  uint32_t a_mar = conf.a_ + conf.ab_ + conf.aB_;
  uint32_t A_mar = conf.A_ + conf.Ab_ + conf.AB_;
  uint32_t b_mar = conf.b_ + conf.ab_ + conf.Ab_;
  uint32_t B_mar = conf.B_ + conf.aB_ + conf.AB_;
  uint32_t degree = a_mar + A_mar + b_mar + B_mar;

  assert(a_mar + A_mar > 0 && b_mar + B_mar > 0);

  return table[degree][a_mar][A_mar][b_mar]
              [conf.ab_][conf.aB_][conf.Ab_][conf.AB_];
}
