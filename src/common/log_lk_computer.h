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

#ifndef LDHELMET_COMMON_LOG_LK_COMPUTER_H_
#define LDHELMET_COMMON_LOG_LK_COMPUTER_H_

#include <string>
#include <vector>

#include "common/change_point.h"
#include "common/site_map_log_lk.h"

struct CompareSnpID {
  inline bool operator()(ChangePoint const &a, ChangePoint const &b) const {
    return a.snp_id_ < b.snp_id_;
  }
};

inline std::string show(std::vector<ChangePoint> const &change_points) {
  std::string change_point_string;
  for (std::vector<ChangePoint>::const_iterator iter = change_points.begin();
       iter != change_points.end();
       ++iter) {
    change_point_string += iter->GetString();
  }
  return change_point_string;
}

class LogLkComputer {
 public:
  explicit LogLkComputer(SiteMapLogLk const &site_map_log_lk, bool interpolate);

  double ComputeLogLkFromSites(double genetic_distance,
                               size_t site0,
                               size_t site1) const;
  inline double GetLogLk(size_t site0, size_t site1, size_t rho_id) const {
    assert(site0 < site1);
    return site_map_log_lk_.site_log_lks_[rho_id * site_map_log_lk_.rho_stride_
                                        + site0 * site_map_log_lk_.window_size_
                                        + (site1 - site0)];
  }
  SiteMapLogLk const &site_map_log_lk_;
 private:
  bool const interpolate_;
};

class FullLogLkComputer {
 public:
  FullLogLkComputer(std::vector<uint64_t> const &snp_pos,
                    LogLkComputer const &log_lk_computer);

  double ComputeLogLk(std::vector<double> const &cum_rho_map) const;

  std::vector<double> ConstructCumRhoMap(
      std::vector<ChangePoint> const &change_points) const;

  double ComputeLogLkConstantRho(double rate) const;

  double ComputeMaxLkConstantRho(double max_lk_start,
                                 double max_lk_end,
                                 double max_lk_resolution) const;

  std::vector<uint64_t> const &snp_pos_;
  LogLkComputer const log_lk_computer_;
};

#endif  // LDHELMET_COMMON_LOG_LK_COMPUTER_H_
