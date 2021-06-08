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

#include "rjmcmc/r_task.h"

#include <stdint.h>
#include <stdio.h>

#include <string>
#include <utility>
#include <vector>

#include <boost/thread.hpp>

#include "common/lk_pade_table.h"
#include "common/mut_mat_prior.h"
#include "rjmcmc/acceptance_log.h"
#include "rjmcmc/handle_output.h"
#include "rjmcmc/priors.h"
#include "rjmcmc/proposals.h"
#include "rjmcmc/ran_num_gen.h"
#include "rjmcmc/rjmcmc.h"

RTaskCounter::RTaskCounter(
    FILE *fp,
    uint64_t max_counter,
    std::vector<std::vector<std::vector<std::pair<uint64_t, double> > > >
        &samples)
    : counter_(),
      max_counter_(max_counter),
      samples_(samples),
      partitions_done(max_counter),
      fp_(fp),
      cur_partition_to_write_() { }

void RTaskCounter::operator()(uint64_t partition_id) {
  boost::lock_guard<boost::mutex> lock(access_);
  assert(cur_partition_to_write_ <= partition_id);

  printf("MCMC done for partition %d (%d/%d).\n",
         static_cast<int>(partition_id),
         static_cast<int>(counter_) + 1,
         static_cast<int>(max_counter_));

  ++counter_;

  assert(partitions_done.at(partition_id) == false);
  partitions_done[partition_id] = true;

  // Write samples_ to output file
  // and check if more partitions should be written out.
  while (cur_partition_to_write_ < partitions_done.size() &&
         partitions_done.at(cur_partition_to_write_)) {
    printf("Writing partition %d to file.\n",
           static_cast<int>(cur_partition_to_write_));

    WriteSamplesToFile(fp_, samples_[cur_partition_to_write_]);

    // Clear vector.
    std::vector<std::vector<std::pair<uint64_t, double> > >().swap(
        samples_[cur_partition_to_write_]);
    ++cur_partition_to_write_;
  }

  fflush(fp_);
}

RTask::RTask(
    RTaskCounter *r_task_counter,
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
    std::pair<uint64_t, uint64_t> const &snp_partition)
    : r_task_counter_(r_task_counter),
      partition_id_(partition_id),
      result_store_(result_store),
      sample_store_(sample_store),
      lk_table_(lk_table),
      snp_seqs_(snp_seqs),
      mut_mat_(mut_mat),
      prior_(prior),
      seed_(seed),
      snp_pos_(snp_pos),
      num_iter_(num_iter),
      burn_in_(burn_in),
      proposal_dist_(proposal_dist),
      acceptance_ratio_(acceptance_ratio),
      max_lk_start_(max_lk_start),
      max_lk_end_(max_lk_end),
      max_lk_resolution_(max_lk_resolution),
      window_size_(window_size),
      stats_thin_(stats_thin),
      snp_partition_overlap_(snp_partition_overlap),
      snp_partition_(snp_partition) { }

void RTask::operator()() {
  assert(snp_partition_overlap_.second
       - snp_partition_overlap_.first + 1 >= 2);

  printf("Preprocessing partition %d.\n",
         static_cast<int>(partition_id_));

  SiteMapLogLk site_map_log_lk(lk_table_,
                               window_size_,
                               snp_seqs_,
                               snp_partition_overlap_.first,
                               snp_partition_overlap_.second,
                               mut_mat_,
                               prior_);

  LogLkComputer log_lk_computer(site_map_log_lk,lk_table_.interpolate_);
  
  // Construct new snp_pos_ to account for snp partitioning.
  std::vector<uint64_t>
    offset_snp_pos(snp_pos_.begin()+snp_partition_overlap_.first,
                   snp_pos_.begin()+snp_partition_overlap_.second);

  assert(snp_partition_.second - snp_partition_.first >= 2);
  assert(snp_partition_.first >= snp_partition_overlap_.first);
  assert(snp_partition_.second >= snp_partition_overlap_.first);

  uint64_t partition_start = snp_partition_.first
                           - snp_partition_overlap_.first;
  uint64_t partition_end = snp_partition_.second
                         - snp_partition_overlap_.first;

  assert(partition_end - partition_start >= 2);
  assert(partition_end - partition_start <=
         snp_partition_overlap_.second - snp_partition_overlap_.first);
  assert(partition_end - partition_start <= offset_snp_pos.size());

  Rjmcmc model(seed_,
               offset_snp_pos,
               log_lk_computer,
               num_iter_,
               burn_in_,
               proposal_dist_,
               acceptance_ratio_,
               max_lk_start_,
               max_lk_end_,
               max_lk_resolution_,
               window_size_,
               stats_thin_,
               partition_start,
               partition_end,
               sample_store_);
  
 printf("MCMC run for partition %d.\n", static_cast<int>(partition_id_));

  model.run();

  // Store result in result_store_.
  result_store_ = model.GetRhoMap();

  (*r_task_counter_)(partition_id_);
}
