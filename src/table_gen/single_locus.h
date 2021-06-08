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

#ifndef LDHELMET_TABLE_GEN_SINGLE_LOCUS_H_
#define LDHELMET_TABLE_GEN_SINGLE_LOCUS_H_

#include <stdint.h>

#include "common/conf.h"

// Compute likelihood at a single locus.
double SingleLocus(double theta, uint32_t a, uint32_t A);

// Take product of two single locus likelihoods.
double ProductSingleLocus(double theta, Conf const &conf);

#endif  // LDHELMET_TABLE_GEN_SINGLE_LOCUS_H_
