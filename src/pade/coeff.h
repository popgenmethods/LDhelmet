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

#ifndef LDHELMET_PADE_COEFF_H_
#define LDHELMET_PADE_COEFF_H_

#include <stdint.h>

#include <vector>

#include "pade/subtable.h"

struct Coeffs {
  std::vector<double> cont_frac;
  std::vector<std::vector<double> > poly_numerator;
  std::vector<std::vector<double> > poly_denominator;
};

// Numerator is indexed by two numbers.
// The first number denotes how many coefficients to "remove" from the Pade
// approximant.
// The second number denotes a particular root given that some number of roots
// have been "removed" from the Pade approximant.
// e.g.,
//   The first root of the numerator given all coefficients are used is:
//     numerator[0][0].
//   The list of roots for the numerator given all coefficients are used is:
//     numerator[0].
//   The list of roots for the numerator given the last coefficient is
//     removed is numerator[1].
//   The list of roots for the numerator given the last two coefficients are
//     removed is numerator[2].
// The same applies to denominator.
struct PadeRoots {
  // [coeff level starting from _last_ coefficient][root_id]
  std::vector<std::vector<double> > numerator;
  // [coeff level starting from _last_ coefficient][root_id]
  std::vector<std::vector<double> > denominator;
};

int ComputeCombinatorialFactor(Conf const &c, Conf const &r);

// Compute q value for a given m and u.
double ComputeQValueIncrementally(Table const &table, int cur_level,
                                  int m, Conf const &c);

// Compute continued fraction coefficients.
// First find coefficients of polynomials using the recursion in
// Baker & Graves-Morris p132.
// The row poly_numerator(k, :) contains the coefficients
// A0, A1,..., Ab of the numerator.
// The continued fraction coefficients are A[0][0], A[1][0], ...
Coeffs ComputeCoeffsForConf(uint32_t num_coeffs,
                            std::vector<double> const &q_values_for_conf,
                            Conf const &conf);

PadeRoots ComputeRelevantRootsForConf(uint64_t num_coeffs,
                                      double defect_threshold,
                                      Coeffs const &coeffs);

#endif  // LDHELMET_PADE_COEFF_H_
