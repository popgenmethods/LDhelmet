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

#include "pade/subtable.h"
#include "pade/subtable-inl.h"

#include <algorithm>

#include "common/conf_gen.h"
#include "common/vector_definitions.h"


double GetTable(Table const &table, uint32_t m, uint32_t u, Conf const &conf) {
  // If m + u = l is odd, then g_u^(m)(a,b,r) = 0.
  if ((m + u) % 2 == 1) {
    return 0.0;
  }
  uint32_t degree = conf.a_ + conf.A_ + conf.b_ + conf.B_
                  + 2*(conf.ab_ + conf.aB_ + conf.Ab_ + conf.AB_);
  assert(degree >= 2*m && degree >= 1);
  uint32_t offset_degree = degree - std::max(static_cast<uint32_t>(1), 2 * m);

  return table[m][u][offset_degree][conf.ab_][conf.aB_][conf.Ab_]
                                   [conf.a_][conf.A_][conf.b_];
}

void SetTable(Table *table, uint32_t m, uint32_t u,
              Conf const &conf, double value) {
  // If m + u = l is odd, then g_u^(m)(a,b,r) = 0,
  //   so the table should never be set if m + u is odd.
  assert((m + u) % 2 == 0);
  uint32_t degree = conf.a_ + conf.A_ + conf.b_ + conf.B_
                  + 2*(conf.ab_ + conf.aB_ + conf.Ab_ + conf.AB_);
  assert(degree >= 2*m && degree >= 1);

  uint32_t offset_degree = degree - std::max(static_cast<uint32_t>(1), 2*m);
  (*table)[m][u][offset_degree][conf.ab_][conf.aB_][conf.Ab_]
                               [conf.a_][conf.A_][conf.b_] = value;
}
