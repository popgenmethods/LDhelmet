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

#ifndef LDHELMET_RJMCMC_RJMCMC_H_
#define LDHELMET_RJMCMC_RJMCMC_H_

#include <stdint.h>

#include <utility>
#include <vector>

#include "common/log_lk_computer.h"
#include "common/seq_process.h"
#include "common/binary_search.h"
#include "rjmcmc/acceptance_log.h"
#include "rjmcmc/post_rho_map.h"
#include "rjmcmc/priors.h"
#include "rjmcmc/proposals.h"
#include "rjmcmc/ran_num_gen.h"

typedef std::vector<double> LogLkMap;

class Rjmcmc {
 public:
  Rjmcmc(RanNumGen::SeedType seed,
         std::vector<uint64_t> const &snp_pos,
         LogLkComputer const &log_lk_computer,
         uint32_t num_iter,
         uint32_t burn_in,
         ProposalDist const &proposal_dist,
         AcceptanceRatio const &acceptance_ratio,
         double max_lk_start,
         double max_lk_end,
         double max_lk_resolution,
         uint32_t window_size,
         uint64_t stats_thin,
         uint64_t partition_start,
         uint64_t partition_end,
         std::vector<std::vector<std::pair<uint64_t, double> > > &sample_store);

  // Run iterations.
  void run();

  // Update MCMC state.
  void Update();

  // Change rate of block.
  void PerformChange();

  // Extend endpoint of block.
  void PerformExtend();

  // Split block into two blocks.
  void PerformSplit();

  // Merge two blocks into one block.
  void PerformMerge();

  // Updates proposed_log_lk_map_.
  double ProposeLogLk(size_t left_snp_id, size_t right_snp_id);

  // Updates the log likelihood map.
  void UpdateLogLkMap(size_t left_snp_id, size_t right_snp_id);


  // Updates cum_rho_map_.
  void UpdateCumRhoMap(std::vector<double> *in_cum_rho_map,
                       size_t begin_change_point_id);

  void CopyCumRhoMap(std::vector<double> *source,
                     std::vector<double> *dest,
                     size_t start_snp_id) const;

  inline std::vector<double> const& GetRhoMap() const {
    return post_rho_map_.GetRhoMap();
  }

  // Records sample from MCMC, handling partitions appropriately.
  void RecordSample();

  RanNumGen rng_;
  RanNumGen::UniformGenType uniform_gen_;  // uniform[0,1] rng.

  // For proposing a change in rate to a block; uniform[-1/2,1/2] rng.
  RanNumGen::UniformGenType rate_change_gen_;

  // Parameters.
  uint32_t const num_iter_;  // Number of iterations to run, excluding burn-in.
  uint32_t const burn_in_;   // Number of itreations for burn-in.

  uint32_t const window_size_;

  uint64_t const stats_thin_;  // Interval in between recording samples.

  AcceptanceRatio const &acceptance_ratio_;

  Proposals proposals_;

  // Data.
  std::vector<uint64_t> const &snp_pos_;
  LogLkComputer const log_lk_computer_;

  // State.
  uint64_t iteration_id_;
  std::vector<ChangePoint> change_points_;
  std::vector<double> cum_rho_map_;
  std::vector<double> proposed_cum_rho_map_;

  LogLkMap log_lk_map_;
  LogLkMap proposed_log_lk_map_;

  double cur_log_lk;

  // Output.
  PostRhoMap post_rho_map_;

  // Debug variables.
  mutable uint32_t rechosen_snp_;

  mutable bool burn_in_p_;
  mutable AcceptanceLog accept_log_burn_in_;
  mutable AcceptanceLog accept_log_run_;

  // Other variables.
  uint64_t partition_start_, partition_end_;

  std::vector<std::vector<std::pair<uint64_t, double> > > &sample_store_;
};

#endif  // LDHELMET_RJMCMC_RJMCMC_H_
