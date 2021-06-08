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

#include "post_to_text/post_to_text.h"

#include <stdio.h>
#include <stdint.h>

#include <algorithm>
#include <string>

#include <boost/lexical_cast.hpp>

#include "common/conf.h"
#include "common/version_number.h"

boost::tuple<RjmcmcParams, RjmcmcResults> ReadPostData(FILE *fp) {
  RjmcmcParams params;
  RjmcmcResults results;

  // Read header.
  int num_read;
  num_read = fread(reinterpret_cast<char *>(&params.version_bit_string_),
                   sizeof(params.version_bit_string_), 1, fp);
  assert(num_read == 1);

  uint64_t version_number = RJMCMC_VERSION_OUTPUT;
  uint64_t version_bit_string = RJMCMC_SALT + version_number;
  if (params.version_bit_string_ != version_bit_string) {
    fprintf(stderr,
            "Input file is not in proper format or has the wrong version "
            "number.\n");
    std::exit(1);
  }

  num_read = fread(reinterpret_cast<char *>(&params.seed_),
                   sizeof(params.seed_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.burn_in_),
                   sizeof(params.burn_in_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.num_iter_),
                   sizeof(params.num_iter_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.stats_thin_),
                   sizeof(params.stats_thin_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.block_penalty_),
                   sizeof(params.block_penalty_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.prior_rate_mean_),
                   sizeof(params.prior_rate_mean_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.window_size_),
                   sizeof(params.window_size_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.partition_length_),
                   sizeof(params.partition_length_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.overlap_length_),
                   sizeof(params.overlap_length_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.pade_resolution_),
                   sizeof(params.pade_resolution_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.pade_max_),
                   sizeof(params.pade_max_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.max_lk_start_),
                   sizeof(params.max_lk_start_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.max_lk_end_),
                   sizeof(params.max_lk_end_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.max_lk_resolution_),
                   sizeof(params.max_lk_resolution_), 1, fp);
  assert(num_read == 1);

  num_read = fread(reinterpret_cast<char *>(&params.num_snps_),
                   sizeof(params.num_snps_), 1, fp);
  assert(num_read == 1);

  for (size_t i = 0; i < params.num_snps_; ++i) {
    uint64_t snp_pos1;
    num_read = fread(reinterpret_cast<char *>(&snp_pos1),
                     sizeof(snp_pos1), 1, fp);
    assert(num_read == 1);
    params.snp_pos_.push_back(snp_pos1);
  }
  assert(params.snp_pos_.size() == params.num_snps_);

  num_read = fread(reinterpret_cast<char*>(&params.num_snp_partitions_),
                   sizeof(params.num_snp_partitions_), 1, fp);
  assert(num_read == 1);

  for (size_t i = 0; i < params.num_snp_partitions_; ++i) {
    uint64_t begin, end;
    num_read = fread(reinterpret_cast<char *>(&begin), sizeof(begin), 1, fp);
    assert(num_read == 1);
    num_read = fread(reinterpret_cast<char *>(&end), sizeof(end), 1, fp);
    assert(num_read == 1);
    params.snp_partitions_.push_back(std::make_pair(begin, end));
  }
  assert(params.snp_partitions_.size() == params.num_snp_partitions_);

  for (size_t i = 0; i < params.num_snp_partitions_; ++i) {
    uint64_t begin, end;
    num_read = fread(reinterpret_cast<char *>(&begin), sizeof(begin), 1, fp);
    assert(num_read == 1);
    num_read = fread(reinterpret_cast<char *>(&end), sizeof(end), 1, fp);
    assert(num_read == 1);
    params.snp_partitions_overlap_.push_back(std::make_pair(begin, end));
  }
  assert(params.snp_partitions_overlap_.size() == params.num_snp_partitions_);

  // Compute number of samples.
  if (params.num_iter_ > 0) {
    params.num_samps_ = 1 + (params.num_iter_ - 1) / params.stats_thin_;
  } else {
    params.num_samps_ = 0;
  }

  // Read posterior mean.
  for (size_t i = 0; i < params.snp_pos_.size(); ++i) {
    double rate;
    num_read = fread(reinterpret_cast<char *>(&rate), sizeof(rate), 1, fp);
    assert(num_read == 1);
    results.mean_.push_back(rate);
  }
  assert(results.mean_.back() == -1.0);

  return boost::make_tuple(params, results);
}

double FindPerc(std::vector<double> const &xs, double perc) {
  if (perc < 0.0 || perc > 1.0) {
    fprintf(stderr,
            "Percentile must be between 0.0 and 1.0, inclusive.\n");
    std::exit(1);
  }

  if (xs.empty()) {
    fprintf(stderr,
            "Array must be non-empty.\n");
    std::exit(1);
  }

  assert(xs.size() > 0);

  std::vector<double> xs2 = xs;
  std::sort(xs2.begin(), xs2.end());
  double p = (xs.size() - 1) * perc;
  size_t n0 = static_cast<size_t>(p);
  assert(n0 < xs2.size());
  assert(static_cast<double>(n0) <= p);

  double ret;
  if (n0 == xs2.size() - 1) {  // Special case when perc == 1.0.
    ret = xs2[n0];
  } else {
    assert(n0 + 1 < xs2.size());
    ret = xs2[n0] * (1.0 - (p - n0)) + xs2[n0 + 1]*(p - n0);
  }

  return ret;
}

std::string Show(RjmcmcParams const &params) {
  std::string tmp;
  tmp += "version_bit_string:";
  tmp += boost::lexical_cast<std::string>(params.version_bit_string_);
  tmp += ',';
  tmp += "seed:";
  tmp += boost::lexical_cast<std::string>(params.seed_);
  tmp += ',';
  tmp += "burn_in:";
  tmp += boost::lexical_cast<std::string>(params.burn_in_);
  tmp += ',';
  tmp += "num_iter:";
  tmp += boost::lexical_cast<std::string>(params.num_iter_);
  tmp += ',';
  tmp += "stats_thin:";
  tmp += boost::lexical_cast<std::string>(params.stats_thin_);
  tmp += ',';
  tmp += "block_penalty:";
  tmp += boost::lexical_cast<std::string>(params.block_penalty_);
  tmp += ',';
  tmp += "prior_rate_mean:";
  tmp += boost::lexical_cast<std::string>(params.prior_rate_mean_);
  tmp += ',';
  tmp += "window_size:";
  tmp += boost::lexical_cast<std::string>(params.window_size_);
  tmp += ',';
  tmp += "partition_length:";
  tmp += boost::lexical_cast<std::string>(params.partition_length_);
  tmp += ',';
  tmp += "overlap_length:";
  tmp += boost::lexical_cast<std::string>(params.overlap_length_);
  tmp += ',';
  tmp += "pade_resolution:";
  tmp += boost::lexical_cast<std::string>(params.pade_resolution_);
  tmp += ',';
  tmp += "pade_max:";
  tmp += boost::lexical_cast<std::string>(params.pade_max_);
  tmp += ',';
  tmp += "max_lk_start:";
  tmp += boost::lexical_cast<std::string>(params.max_lk_start_);
  tmp += ',';
  tmp += "max_lk_resolution:";
  tmp += boost::lexical_cast<std::string>(params.max_lk_resolution_);
  tmp += ',';
  tmp += "num_snps:";
  tmp += boost::lexical_cast<std::string>(params.num_snps_);
  tmp += ',';
  tmp += "num_snp_partitions:";
  tmp += boost::lexical_cast<std::string>(params.num_snp_partitions_);
  tmp += ',';
  tmp += "num_samps:";
  tmp += boost::lexical_cast<std::string>(params.num_samps_);

  return tmp;
}
