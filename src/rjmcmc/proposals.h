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

#ifndef LDHELMET_RJMCMC_PROPOSALS_H
#define LDHELMET_RJMCMC_PROPOSALS_H

#include <string>

struct ProposalDist {
  double change;  // Change rate.
  double extend;  // Extend block.
  double split;   // Split block.
  double merge;   // Merge blocks.

  ProposalDist() { }

  ProposalDist(double change, double extend, double split, double merge)
      : change(change), extend(extend), split(split), merge(merge) { }
};

std::string Show(ProposalDist const &prop);

class Proposals {
 public:
  ProposalDist dist_;
  ProposalDist cum_;
  explicit Proposals(ProposalDist const &proposal_dist);
};

#endif  // LDHELMET_RJMCMC_PROPOSALS_H
