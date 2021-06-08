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

#ifndef LDHELMET_RJMCMC_POST_RHO_MAP_H_
#define LDHELMET_RJMCMC_POST_RHO_MAP_H_

#include <stddef.h>
#include <stdint.h>

#include <vector>

class PostRhoMap {
 public:
  explicit PostRhoMap(size_t num_snps);

  void Update(std::vector<uint64_t> const &snp_pos,
              std::vector<double> const &cum_rho_map);

  std::vector<double> ComputeRhoMap(
      std::vector<uint64_t> const &snp_pos,
      std::vector<double> const &cum_rho_map) const;

  // Compute the cumulative rho map from the posterior mean map.
  std::vector<double> GetCumRhoMap(
      std::vector<uint64_t> const &snp_pos_) const;

  double BackgroundRho(std::vector<uint64_t> const &snp_pos_) const;

  inline std::vector<double> const &GetRhoMap() const {
    return map_;
  }

  uint64_t num_iter_;
  std::vector<double> map_;
};

#endif  // LDHELMET_RJMCMC_POST_RHO_MAP_H_
