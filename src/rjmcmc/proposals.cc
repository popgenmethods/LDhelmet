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

#include "rjmcmc/proposals.h"

#include <stdio.h>

#include <cassert>
#include <string>

#include <boost/lexical_cast.hpp>

std::string Show(ProposalDist const &prop) {
  return std::string("[")
       + "change: " + boost::lexical_cast<std::string>(prop.change) + ", "
       + "extend: " + boost::lexical_cast<std::string>(prop.extend) + ", "
       + "split: " + boost::lexical_cast<std::string>(prop.split) + ", "
       + "merge: " + boost::lexical_cast<std::string>(prop.merge)
       + "]";
}

Proposals::Proposals(ProposalDist const &proposal_dist) : dist_(proposal_dist) {
  // Proposal probabilities must be non-negative.
  if (!(dist_.change >= 0.0 &&
        dist_.extend >= 0.0 &&
        dist_.split  >= 0.0 &&
        dist_.merge  >= 0.0)) {
    fprintf(stderr,
            "Proposal distribution has negative components.\n");
    std::exit(1);
  }

  // Proposal probabilities must have at least one positive entry.
  if (dist_.change <= 0.0 &&
      dist_.extend <= 0.0 &&
      dist_.split  <= 0.0 &&
      dist_.merge  <= 0.0) {
    fprintf(stderr,
            "Proposal distribution has no positive components.\n");
    std::exit(1);
  }

  if ((dist_.split > 0.0 && dist_.merge == 0.0) ||
      (dist_.merge > 0.0 && dist_.split == 0.0)) {
    fprintf(stderr,
      "Transition probabilities are not reversible due to improper "
      "split and merge proposal probabilities. "
      "They must either be both zero or both positive.\n");
    std::exit(1);
  }

  // Normalize proposal distribution.
  double prop_norm = dist_.change + dist_.extend + dist_.split + dist_.merge;

  assert(prop_norm > 0.0);

  dist_.change /= prop_norm;
  dist_.extend /= prop_norm;
  dist_.split /= prop_norm;
  dist_.merge /= prop_norm;

  // Compute cumulative proposal distribution.
  cum_.change = dist_.change;
  cum_.extend = cum_.change + dist_.extend;
  cum_.split = cum_.extend + dist_.split;
  cum_.merge = 1.0;  // cum_.split + dist_.merge should be 1.0.
}
