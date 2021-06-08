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

#include <stdint.h>
#include <stdio.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include "common/lk_pade_table.h"
#include "common/load_data.h"
#include "common/log_lk_computer.h"
#include "common/mut_mat_prior.h"
#include "common/site_map_log_lk.h"
#include "common/version_number.h"
#include "max_lk/max_lk_options.h"

int MaxLkComponent(std::string const base_command, int argc, char **argv) {
  std::string version_string = MAXLK_VERSION_STRING;
  CmdLineOptionsMaxLk cmd_line_options(base_command,
                                       argc, argv, version_string);
  if (!cmd_line_options.success()) {
    std::exit(1);
  }

  printf("%s\n", ShowCommandLineOptions(argc, argv).c_str());

  printf("Number of threads: %d.\n",
         static_cast<int>(cmd_line_options.num_threads_));
  assert(cmd_line_options.num_threads_ > 0);

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

  printf("Loading likelihood and Pade files.\n");

  LkTable lk_table =
    LoadLikelihoodAndPade(cmd_line_options.lk_file_,
                          cmd_line_options.pade_file_,
                          cmd_line_options.pade_resolution_,
                          cmd_line_options.pade_max_,
                          cmd_line_options.window_size_,
                          snp_seqs,
                          cmd_line_options.num_threads_);

  printf("Preprocessing sequence data.\n");

  size_t snp_begin = 0, snp_end = snp_seqs.front().size();
  SiteMapLogLk site_map_log_lk(lk_table,
                               cmd_line_options.window_size_,
                               snp_seqs,
                               snp_begin,
                               snp_end,
                               mut_mat,
                               prior);

  LogLkComputer log_lk_computer(site_map_log_lk, lk_table.interpolate_);
  FullLogLkComputer full_log_lk_computer(snp_pos, log_lk_computer);

  printf("Computing maximum likelihood estimate.\n");

  double max_lk_rate =
    full_log_lk_computer.ComputeMaxLkConstantRho(
        cmd_line_options.max_lk_start_,
        cmd_line_options.max_lk_end_,
        cmd_line_options.max_lk_resolution_);

  printf("Maximum likelihood estimate: %.8f.\n", max_lk_rate);

  return 0;
}
