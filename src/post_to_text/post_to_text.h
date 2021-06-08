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

#ifndef LDHELMET_POST_TO_TEXT_POST_TO_TEXT_H_
#define LDHELMET_POST_TO_TEXT_POST_TO_TEXT_H_

#include <stdint.h>
#include <stdio.h>

#include <string>
#include <utility>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include "common/conf.h"

class RjmcmcParams {
 public:
  uint64_t version_bit_string_;
  uint64_t seed_;

  uint64_t burn_in_;
  uint64_t num_iter_;

  uint64_t stats_thin_;
  double block_penalty_;
  double prior_rate_mean_;
  uint64_t window_size_;

  uint64_t partition_length_;
  uint64_t overlap_length_;

  double pade_resolution_;
  double pade_max_;

  double max_lk_start_;
  double max_lk_end_;
  double max_lk_resolution_;

  uint64_t num_snps_;
  std::vector<uint64_t> snp_pos_;

  uint64_t num_snp_partitions_;
  std::vector<std::pair<uint64_t, uint64_t> > snp_partitions_;
  std::vector<std::pair<uint64_t, uint64_t> > snp_partitions_overlap_;

  uint64_t num_samps_;
};

class RjmcmcResults {
 public:
  std::vector<double> mean_;
};

boost::tuple<RjmcmcParams, RjmcmcResults> ReadPostData(FILE *fp);

double FindPerc(std::vector<double> const &xs, double perc);

std::string Show(RjmcmcParams const &params);

#endif  // LDHELMET_POST_TO_TEXT_POST_TO_TEXT_H_
