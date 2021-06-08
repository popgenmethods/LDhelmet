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

#ifndef LDHELMET_MAX_LK_MAX_LK_OPTIONS_H_
#define LDHELMET_MAX_LK_MAX_LK_OPTIONS_H_

// Header file implementing a basic command line and configuration file parser.

#include <string>

#include "common/command_line_options.h"

class CmdLineOptionsMaxLk : public CmdLineOptions {
 public:
  CmdLineOptionsMaxLk(std::string const base_command,
                      int argc,
                      char **argv,
                      std::string version);

  uint32_t num_threads_;

  uint32_t window_size_;      // Window size.
  std::string seq_file_;      // Snp sequence file.
  std::string lk_file_;       // Likelihood table file.
  std::string pade_file_;     // Pade coefficient file.
  std::string prior_file_;    // Prior on ancestral allele.
  std::string mut_mat_file_;  // Mutation matrix.

  double pade_resolution_;    // Pade grid increment.
  double pade_max_;           // Max value in Pade grid.

  std::string pos_file_;      // Snp positions file for alternative input
                              //   format.
  std::string snps_file_;     // Snps file for alternative input format.

  double max_lk_start_;       // Maximum likelihood estimate start.
  double max_lk_end_;         // Maximum likelihood estimate end.
  double max_lk_resolution_;  // Maximum likelihood estimate increment.

 private:
  std::string GetUsageString(std::string const base_command,
                           std::string const program_name) const;

  bool ProperInput(
    boost::program_options::variables_map const &variables_map) const;
};

#endif  // LDHELMET_MAX_LK_MAX_LK_OPTIONS_H_
