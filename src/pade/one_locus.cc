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

#include "pade/one_locus.h"

#include <stddef.h>
#include <stdio.h>

#include <cstdlib>

double OneLocusSF(double theta, std::vector<int> const &n) {
  int nsum = 0;
  for (size_t i = 0; i < n.size(); ++i) {
    nsum += n[i];
  }

  double prob = 0.0;

  if (n[0] != nsum) {
    for (int la = 1; la <= n[0]; ++la) {
      double temp = 1.0;
      for (int i = 1; i < la; ++i) {
        temp *= n[0] - i;
      }
      for (int i = 0; i < nsum - 1 - la; ++i) {
        temp *= nsum - 1 - la - i;
      }
      temp *= la / (la + theta);
      prob += temp;
    }

    prob *= theta;
    for (int i = 0; i < nsum - 1; ++i) {
      prob /= 1 + theta + i;
    }
  } else {
    prob = 1.0;
    for (int i = 0; i < nsum - 1; ++i) {
      prob *= (nsum - 1 - i) / (1 + theta + i);
    }
  }

  if (prob > 1.0) {
    fprintf(stderr, "Error: Probability greater than 1.\n");
    std::exit(1);
  }

  for (int i = 2; i <= nsum; ++i) {
    prob /= i;
    for (size_t j = 0; j < n.size(); ++j) {
      if (i <= n[j]) {
        prob *= i;
      }
    }
  }

  return prob;
}
