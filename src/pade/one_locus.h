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

#ifndef LDHELMET_PADE_ONE_LOCUS_H_
#define LDHELMET_PADE_ONE_LOCUS_H_

#include <vector>

// This function computes the one-locus sampling probability jointly with at
// most one mutation event at this locus, assuming that allele 1 is ancestral.
double OneLocusSF(double theta, std::vector<int> const &n);

#endif  // LDHELMET_PADE_ONE_LOCUS_H_
