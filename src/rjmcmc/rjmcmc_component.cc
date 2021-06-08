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

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include <boost/asio.hpp>
#include <boost/lexical_cast.hpp>

#include "common/load_data.h"
#include "common/seq_file_parse.h"
#include "common/seq_process.h"
#include "common/version_number.h"
#include "rjmcmc/handle_output.h"
#include "rjmcmc/r_task.h"
#include "rjmcmc/ran_num_gen.h"
#include "rjmcmc/rjmcmc.h"
#include "rjmcmc/rjmcmc_options.h"

int RjmcmcComponent(std::string const base_command, int argc, char **argv) {
  uint64_t version_number = RJMCMC_VERSION_OUTPUT;
  uint64_t version_bit_string = RJMCMC_SALT + version_number;

  std::string version_string = RJMCMC_VERSION_STRING;
  CmdLineOptionsRjmcmc cmd_line_options(base_command,
                                        argc,
                                        argv,
                                        version_string);

  if (!cmd_line_options.success()) {
    std::exit(1);
  }

  printf("%s\n", ShowCommandLineOptions(argc, argv).c_str());

  printf("Number of threads: %d.\n",
         static_cast<int>(cmd_line_options.num_threads_));
  assert(cmd_line_options.num_threads_ > 0);

  printf("Window size: %d.\n",
         static_cast<int>(cmd_line_options.window_size_));

  MutationMatrix mut_mat;
  std::vector<std::string> snp_seqs;
  std::vector<uint64_t> snp_pos;
  Prior prior;
  boost::tie(mut_mat, snp_seqs, snp_pos, prior)
    = LoadData(cmd_line_options.num_threads_,
               cmd_line_options.mut_mat_file_,
               cmd_line_options.seq_file_,
               cmd_line_options.pos_file_,
               cmd_line_options.snps_file_,
               cmd_line_options.prior_file_);

  // Compute seq partitions.
  SnpPartitions snp_partitions;
  SnpPartitions snp_partitions_overlap;
  boost::tie(snp_partitions, snp_partitions_overlap) =
      ComputeSnpPartitions(snp_pos,
                           cmd_line_options.partition_length_,
                           cmd_line_options.overlap_length_);

  assert(snp_partitions.size() > 0);
  assert(snp_partitions.size() == snp_partitions_overlap.size());

  printf("Number of partitions: %d.\n",
         static_cast<int>(snp_partitions.size()));

  printf("Partitions: %s.\n", ShowSnpPartitions(snp_partitions).c_str());

  FILE *fp = fopen(cmd_line_options.output_file_.c_str(), "wb");

  int64_t marg_exp_pos = WriteHeaderToFile(fp,
                                           version_bit_string,
                                           cmd_line_options,
                                           snp_pos,
                                           snp_partitions,
                                           snp_partitions_overlap);
  
  LkTable lk_table = LoadLikelihoodAndPade(cmd_line_options.lk_file_,
                                           cmd_line_options.pade_file_,
                                           cmd_line_options.pade_resolution_,
                                           cmd_line_options.pade_max_,
                                           cmd_line_options.window_size_,
                                           snp_seqs,
                                           cmd_line_options.num_threads_);

  uint32_t const num_iter = cmd_line_options.num_iter_;
  uint32_t const burn_in = cmd_line_options.burn_in_;
  double const block_penalty = cmd_line_options.block_penalty_;
  double const prior_rate_mean = cmd_line_options.prior_rate_mean_;

  // Proposals: change, extend, split, merge.
  ProposalDist proposal_dist(0.4, 0.2, 0.2, 0.2);

  AcceptanceRatio const &acceptance_ratio
      = ExponentialPrior(block_penalty, prior_rate_mean);

  // Declare variables to store results.
  std::vector<std::vector<double> > results(snp_partitions_overlap.size());
  std::vector<std::vector<std::vector<std::pair<uint64_t, double> > > >
      samples(snp_partitions_overlap.size());

  // Use thread pool.
  {
    RTaskCounter r_task_counter(fp, snp_partitions_overlap.size(), samples);
    RanNumGen seed_generator(cmd_line_options.seed_);
    boost::thread_group threads;
    boost::asio::io_service io_service;

    {
      boost::asio::io_service::work work(io_service);
      for (size_t i = 0; i < cmd_line_options.num_threads_; ++i) {
        threads.create_thread(
            boost::bind(&boost::asio::io_service::run, &io_service));
      }

      for (size_t partition_id = 0;
           partition_id < snp_partitions_overlap.size();
           ++partition_id) {
        RanNumGen::SeedType thread_seed =
          seed_generator.GetRandomNumberGen()();

        std::vector<double> &result_store = results[partition_id];
        std::vector<std::vector<std::pair<uint64_t, double> > >
          &sample_store = samples[partition_id];

        RTask model_task(&r_task_counter,
                         partition_id,
                         result_store,
                         sample_store,
                         lk_table,
                         snp_seqs,
                         mut_mat,
                         prior,
                         thread_seed,
                         snp_pos,
                         num_iter,
                         burn_in,
                         proposal_dist,
                         acceptance_ratio,
                         cmd_line_options.max_lk_start_,
                         cmd_line_options.max_lk_end_,
                         cmd_line_options.max_lk_resolution_,
                         cmd_line_options.window_size_,
                         cmd_line_options.stats_thin_,
                         snp_partitions_overlap[partition_id],
                         snp_partitions[partition_id]);

        io_service.post(model_task);
      }
    }

    threads.join_all();
  }

  printf("MCMC complete. Merging results.\n");

  std::vector<double> rho_map = MergeResults(snp_pos.size(),
                                             snp_partitions,
                                             snp_partitions_overlap,
                                             results);

  {
    int ret = fseek(fp, marg_exp_pos, SEEK_SET);
    if (ret) {
      fprintf(stderr, "Error seeking in output file.\n");
      std::exit(1);
    }

    printf("Writing results to %s.\n",
           cmd_line_options.output_file_.c_str());

    WriteResultToFile(fp, snp_pos, rho_map);
  }

  return 0;
}
