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

#include "common/log_lk_computer.h"

#include <stdio.h>

#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

#include "common/change_point.h"
#include "common/site_map_log_lk.h"

LogLkComputer::LogLkComputer(SiteMapLogLk const &site_map_log_lk, bool const interpolate)
    : site_map_log_lk_(site_map_log_lk), interpolate_(interpolate) { }

double LogLkComputer::ComputeLogLkFromSites(double genetic_distance,
                                            size_t site0,
                                            size_t site1) const {
  if (genetic_distance < site_map_log_lk_.rho_finder_.rho_list_[0]) {
    fprintf(stderr,
            "Error: Encountered a genetic distance less than all "
            "rho values in the likelihood table. Either a genetic "
            "distance is negative, which is a bug in the program, or "
            "the likelihood table doesn't include rho value zero.\n");
    std::exit(1);
  }

  size_t left_rho_id = site_map_log_lk_.rho_finder_.GetRhoID(genetic_distance);
  size_t right_rho_id = left_rho_id + 1;
  if(!this->interpolate_){
    if (right_rho_id < site_map_log_lk_.rho_finder_.rho_list_.size()) {
      return GetLogLk(site0, site1, right_rho_id);
    }
    else{
      return GetLogLk(site0, site1, left_rho_id);
    }
  }
  double left_conf_log_lk =  GetLogLk(site0, site1, left_rho_id);
  if(std::isinf(left_conf_log_lk)){
      fprintf(stderr,
              "Encountered a zero likelihood configuration. "
              "This usually indicates an issue with the interpolation\n");
      std::exit(1); 
  }
  double right_conf_log_lk;
  double left_distance = genetic_distance - site_map_log_lk_.rho_finder_.rho_list_[left_rho_id];
  double right_distance;
  if (right_rho_id < site_map_log_lk_.rho_finder_.rho_list_.size()) {
    right_conf_log_lk = GetLogLk(site0, site1, right_rho_id);
    right_distance = site_map_log_lk_.rho_finder_.rho_list_[right_rho_id] - genetic_distance;
  } else {
    assert(right_rho_id == site_map_log_lk_.rho_finder_.rho_list_.size());
    return GetLogLk(site0, site1, right_rho_id - 1);
    
  }
  assert(left_distance + right_distance > 0);
  //interpolate
  double conf_log_lk = (right_distance * left_conf_log_lk + left_distance * right_conf_log_lk)
                                                / (left_distance + right_distance);
  assert(conf_log_lk <= 0.0);
  assert(!std::isinf(conf_log_lk));
  assert(!std::isnan(conf_log_lk));

  return conf_log_lk;
}

FullLogLkComputer::FullLogLkComputer(std::vector<uint64_t> const &snp_pos,
                                     LogLkComputer const &log_lk_computer)
    : snp_pos_(snp_pos), log_lk_computer_(log_lk_computer) {
  assert(snp_pos_.size() == log_lk_computer_.site_map_log_lk_.len_snp_seq_);
}

double FullLogLkComputer::ComputeLogLk(
    std::vector<double> const &cum_rho_map) const {
  size_t begin_snp_id = 0;
  size_t end_snp_id = snp_pos_.size();

  double log_lk = 0.0;
  for (size_t site0 = begin_snp_id; site0 < end_snp_id; ++site0) {
    size_t max_site1 =
      std::min(site0 + log_lk_computer_.site_map_log_lk_.window_size_,
               snp_pos_.size());
    for (size_t site1 = site0 + 1; site1 < max_site1; ++site1) {
      double genetic_distance = cum_rho_map[site1] - cum_rho_map[site0];
      double conf_log_lk =
        log_lk_computer_.ComputeLogLkFromSites(genetic_distance,
                                               site0, site1);
      log_lk += conf_log_lk;
    }
  }
  return log_lk;
}

std::vector<double> FullLogLkComputer::ConstructCumRhoMap(
    std::vector<ChangePoint> const &change_points) const {
  std::vector<double> cum_rho_map(snp_pos_.size(), 0.0);

  double rho_sum = 0.0;
  for (size_t change_point_id = 0;
       change_point_id < change_points.size() - 1;
       ++change_point_id) {
    for (size_t snp_id = change_points[change_point_id].snp_id_;
         snp_id < change_points[change_point_id + 1].snp_id_;
         ++snp_id) {
      assert(snp_id < snp_pos_.size());
      cum_rho_map[snp_id] = rho_sum;
      rho_sum += static_cast<double>(snp_pos_[snp_id + 1] - snp_pos_[snp_id])
               * change_points[change_point_id].rate_;
      assert(change_points[change_point_id].rate_ >= 0.0);
    }
  }
  cum_rho_map.back() = rho_sum;

  return cum_rho_map;
}

double FullLogLkComputer::ComputeLogLkConstantRho(double rate) const {
  std::vector<ChangePoint> change_points;
  change_points.push_back(ChangePoint(0, rate));
  change_points.push_back(ChangePoint(snp_pos_.size(), -1.0));
  double log_lk = ComputeLogLk(ConstructCumRhoMap(change_points));
  return log_lk;
}

double FullLogLkComputer::ComputeMaxLkConstantRho(
    double max_lk_start,
    double max_lk_end,
    double max_lk_resolution) const {
  if (max_lk_start < 0.0 || max_lk_end < 0.0 || max_lk_resolution <= 0.0) {
    fprintf(stderr,
            "Error: Max likelihood start and max likelihood end must both be "
            "non-negative. "
            "Max likelihood resolution must be greater than 0.0.\n");
    std::exit(1);
  }
  if (max_lk_start > max_lk_end) {
    fprintf(stderr,
            "Error: Max likelihood start must be less than or equal to "
            "max likelihood end.\n");
    std::exit(1);
  }
  double max_log_lk = -std::numeric_limits<double>::infinity();
  double max_lk_rate = -1.0;
  uint64_t point_id = 0;
  while (max_lk_start + point_id * max_lk_resolution < max_lk_end) {
    double rate = max_lk_start + point_id*max_lk_resolution;
    double log_lk = ComputeLogLkConstantRho(rate);
    if (log_lk > max_log_lk) {
      max_log_lk = log_lk;
      max_lk_rate = rate;
    }
    ++point_id;
  }
  return max_lk_rate;
}
