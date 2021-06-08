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

#include "common/ncr.h"

#include <algorithm>

uint64_t Ncr(uint64_t n, uint64_t r) {
  if (r > n) {
    return 0;
  }
  uint64_t product = 1;
  for (uint64_t i = 0; i < std::min(r, n - r); ++i) {
    product = product * (n - i) / (i + 1);
  }
  return product;
}
