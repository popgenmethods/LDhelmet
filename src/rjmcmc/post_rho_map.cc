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

#include "rjmcmc/post_rho_map.h"

#include <stddef.h>
#include <stdint.h>

#include <cassert>
#include <vector>

PostRhoMap::PostRhoMap(size_t num_snps) : num_iter_(0),
                                          map_(num_snps - 1, 0.0) { }

void PostRhoMap::Update(std::vector<uint64_t> const &snp_pos,
                        std::vector<double> const &cum_rho_map) {
  assert(map_.size() == cum_rho_map.size() - 1);
  std::vector<double> rho_map = ComputeRhoMap(snp_pos, cum_rho_map);
  for (size_t snp_id = 0; snp_id < map_.size(); ++snp_id) {
    map_[snp_id] =
      (static_cast<double>(num_iter_) * map_[snp_id] + rho_map[snp_id])
      / static_cast<double>(num_iter_ + 1);
  }

  ++num_iter_;
}

std::vector<double> PostRhoMap::ComputeRhoMap(
    std::vector<uint64_t> const &snp_pos,
    std::vector<double> const &cum_rho_map) const {
  std::vector<double> rho_map(cum_rho_map.size() - 1);
  for (size_t snp_id = 0; snp_id < cum_rho_map.size() - 1; ++snp_id) {
    double rate = (cum_rho_map[snp_id + 1] - cum_rho_map[snp_id])
                  / static_cast<double>(snp_pos[snp_id + 1]
                                          - snp_pos[snp_id]);
    rho_map[snp_id] = rate;
  }
  return rho_map;
}

std::vector<double> PostRhoMap::GetCumRhoMap(
    std::vector<uint64_t> const &snp_pos) const {
  assert(snp_pos.size() - 1 == map_.size());
  std::vector<double> cum_rho_map(map_.size() + 1);
  double rho_sum = 0.0;
  for (size_t snp_id = 0; snp_id < snp_pos.size() - 1; ++snp_id) {
    cum_rho_map[snp_id] = rho_sum;
    rho_sum += map_[snp_id]
             * static_cast<double>(snp_pos[snp_id + 1] - snp_pos[snp_id]);
  }
  cum_rho_map.back() = rho_sum;
  return cum_rho_map;
}

double PostRhoMap::BackgroundRho(std::vector<uint64_t> const &snp_pos) const {
  double sum = 0.0;
  for (size_t snp_id = 0; snp_id < map_.size(); ++snp_id) {
    sum += map_[snp_id]
         * static_cast<double>(snp_pos[snp_id + 1] - snp_pos[snp_id]);
  }
  return sum / static_cast<double>(snp_pos.back());
}
