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

#ifndef LDHELMET_RJMCMC_RJMCMC_OPTIONS_H_
#define LDHELMET_RJMCMC_RJMCMC_OPTIONS_H_

// Header file implementing a basic command line and configuration file parser.

#include <string>

#include "common/command_line_options.h"
#include "rjmcmc/ran_num_gen.h"

class CmdLineOptionsRjmcmc : public CmdLineOptions {
 public:
  CmdLineOptionsRjmcmc(std::string const base_command,
                       int argc,
                       char **argv,
                       std::string const &version);
  uint32_t num_threads_;
  RanNumGen::SeedType seed_;

  std::string output_file_;
  uint64_t stats_thin_;     // Interval between recording samples.

  uint32_t num_iter_;       // Number of iterations for RJMCMC.
  uint32_t burn_in_;        // Number of iterations for burn-in
                            //   in addition to num_iter_.

  uint32_t window_size_;    // Window size.
  double block_penalty_;    // Block penalty.
  double prior_rate_mean_;  // Prior rate mean for rates per base.

  std::string seq_file_;      // Snp sequence file.
  std::string lk_file_;       // Likelihood table file.
  std::string pade_file_;     // Pade coefficient file.
  std::string prior_file_;    // Prior on ancestral allele.
  std::string mut_mat_file_;  // Mutation matrix.

  double pade_resolution_;  // Pade grid increment.
  double pade_max_;         // Max value in Pade grid.

  uint64_t partition_length_;  // Length of partitions to break sequence into
  uint64_t overlap_length_;    // amount the partitions overlap each other.
                               // The overlap is added to the partition length.

  std::string pos_file_;    // snp positions file for alternative input
                            //     format
  std::string snps_file_;   // snps file for alternative input format

  double max_lk_start_;
  double max_lk_end_;
  double max_lk_resolution_;

 private:
  std::string GetUsageString(std::string const base_command,
                             std::string const program_name) const;

  bool ProperInput(
      boost::program_options::variables_map const &variablesMap) const;
};

#endif  // LDHELMET_RJMCMC_RJMCMC_OPTIONS_H_
