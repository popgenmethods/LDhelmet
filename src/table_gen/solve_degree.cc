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

#include "table_gen/solve_degree.h"

#include <algorithm>

#include "common/conf.h"
#include "common/mar.h"
#include "common/predicates.h"
#include "common/vector_definitions.h"
#include "table_gen/single_locus.h"
#include "table_gen/table_management.h"

void AssignBaseCases(Vec8 *table,
                     double theta,
                     double rho,
                     uint32_t a_mar,
                     uint32_t A_mar,
                     uint32_t b_mar,
                     uint32_t B_mar) {
  assert((a_mar == 1 && A_mar == 0) || (b_mar == 1 && B_mar == 0));
  Conf conf;
  uint32_t min_ab = std::min(a_mar, b_mar);
  for (conf.ab_ = 0; conf.ab_ <= min_ab; ++conf.ab_) {
    uint32_t min_aB = std::min(a_mar - conf.ab_, B_mar);
    for (conf.aB_= 0; conf.aB_ <= min_aB; ++conf.aB_) {
      uint32_t min_Ab = std::min(b_mar - conf.ab_, A_mar);
      for (conf.Ab_ = 0; conf.Ab_ <= min_Ab; ++conf.Ab_) {
        uint32_t min_AB = std::min(A_mar - conf.Ab_, B_mar - conf.aB_);
        for (conf.AB_ = 0; conf.AB_ <= min_AB; ++conf.AB_) {
          conf.a_ = a_mar - conf.ab_ - conf.aB_;
          conf.A_ = A_mar - conf.Ab_ - conf.AB_;
          conf.b_ = b_mar - conf.ab_ - conf.Ab_;
          conf.B_ = B_mar - conf.aB_ - conf.AB_;

          double prob = ProductSingleLocus(theta, conf);
          SetTable(table, conf, prob);
        }
      }
    }
  }
}
