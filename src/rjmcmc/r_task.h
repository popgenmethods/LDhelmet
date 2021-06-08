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

#ifndef LDHELMET_RJMCMC_R_TASK_H_
#define LDHELMET_RJMCMC_R_TASK_H_

#include <stdint.h>

#include <string>
#include <utility>
#include <vector>

#include <boost/thread.hpp>

#include "common/lk_pade_table.h"
#include "common/mut_mat_prior.h"
#include "rjmcmc/acceptance_log.h"
#include "rjmcmc/priors.h"
#include "rjmcmc/proposals.h"
#include "rjmcmc/ran_num_gen.h"

class RTaskCounter {
 public:
  RTaskCounter(
      FILE *fp,
      uint64_t max_counter,
      std::vector<std::vector<std::vector<std::pair<uint64_t, double> > > >
        &samples);

  void operator()(uint64_t partition_id);

  boost::mutex access_;
  uint64_t counter_;
  uint64_t const max_counter_;
  std::vector<std::vector<std::vector<std::pair<uint64_t, double> > > >
    &samples_;
  std::vector<bool> partitions_done;
  FILE *fp_;
  uint64_t cur_partition_to_write_;
};

class RTask {
 public:
  RTask(RTaskCounter *r_task_counter,
         uint64_t partition_id,
         std::vector<double> &result_store,
         std::vector<std::vector<std::pair<uint64_t, double> > > &sample_store,
         LkTable const &lk_table,
         std::vector<std::string> const &snp_seqs,
         MutationMatrix const &mut_mat,
         Prior const &prior,
         RanNumGen::SeedType seed,
         std::vector<uint64_t> const &snp_pos,
         uint32_t num_iter,
         uint32_t burn_in,
         ProposalDist const &proposal_dist,
         AcceptanceRatio const &acceptance_ratio,
         double max_lk_start,
         double max_lk_end,
         double max_lk_resolution,
         uint32_t window_size,
         uint64_t stats_thin,
         std::pair<uint64_t, uint64_t> const &snp_partition_overlap,
         std::pair<uint64_t, uint64_t> const &snp_partition);

  void operator()();

  RTaskCounter *r_task_counter_;
  uint64_t partition_id_;

  std::vector<double> &result_store_;
  std::vector<std::vector<std::pair<uint64_t, double> > > &sample_store_;

  LkTable const &lk_table_;
  std::vector<std::string> const &snp_seqs_;
  MutationMatrix const &mut_mat_;
  Prior const &prior_;

  RanNumGen::SeedType seed_;
  std::vector<uint64_t> const &snp_pos_;
  uint32_t num_iter_;
  uint32_t burn_in_;
  double block_penalty_;
  ProposalDist const &proposal_dist_;
  AcceptanceRatio const &acceptance_ratio_;

  double max_lk_start_, max_lk_end_, max_lk_resolution_;

  uint32_t window_size_;

  uint64_t stats_thin_;

  std::pair<uint64_t, uint64_t> snp_partition_overlap_;
  std::pair<uint64_t, uint64_t> snp_partition_;
};

#endif  // LDHELMET_RJMCMC_R_TASK_H_
