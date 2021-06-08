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

#include "table_gen/single_locus.h"

#include <stdint.h>

#include "common/conf.h"
#include "common/ncr.h"

double SingleLocus(double theta, uint32_t a, uint32_t A) {
  uint32_t n = a + A;

  assert(n > 0);

  if (n == A) {  // No ancestral allele.
    return 0.0;
  }

  uint32_t combinatorial_factor = Ncr(n, a);

  if (n == a) {  // a - 1 coalescences.
    double prob = 1.0;
    for (uint32_t i = 0; i < n - 1; ++i) {
      prob = prob
             * (static_cast<double>(a) - 1.0 - static_cast<double>(i))
             / (static_cast<double>(n - 1) + theta - static_cast<double>(i));
    }
    return prob / static_cast<double>(combinatorial_factor);
  } else {  // One mutation occurs along with coalescences.
    // Compute ratio of factorial over rising factorial.
    double fac_ratio = 1.0;
    for (uint32_t i = 0; i < n - 1; ++i) {
      fac_ratio = fac_ratio
                 * (static_cast<double>(n - 1) - static_cast<double>(i))
                 / (static_cast<double>(n - 1) - static_cast<double>(i)
                      + theta);
    }

    // Compute sum.
    double sum = 0.0;
    for (uint32_t i = 1; i <= a; ++i) {
      sum += (static_cast<double>(Ncr(a - 1, i - 1))
             / static_cast<double>(Ncr(n - 1, i)))
             * (1.0 / (static_cast<double>(i) + theta));
    }

    return theta * fac_ratio * sum / static_cast<double>(combinatorial_factor);
  }
}

double ProductSingleLocus(double theta, Conf const &conf) {
  return SingleLocus(theta, conf.ComputeaMar(), conf.ComputeAMar())
       * SingleLocus(theta, conf.ComputebMar(), conf.ComputeBMar());
}
