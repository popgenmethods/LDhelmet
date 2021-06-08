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

#ifndef LDHELMET_PADE_COMPUTE_G_H_
#define LDHELMET_PADE_COMPUTE_G_H_

#include <stdint.h>

#include "common/conf.h"
#include "pade/one_locus.h"
#include "pade/subtable.h"

template<class T>
inline double delta(T a, T b) {
  return a == b ? 1.0 : 0.0;
}

double ComputeG(double theta, Table *table,
                uint32_t m, int ind, Conf const &conf);

#endif  // LDHELMET_PADE_COMPUTE_G_H_
