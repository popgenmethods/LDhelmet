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

#include "max_lk/max_lk_options.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <string>

#include "common/command_line_options.h"

static uint32_t const kDefaultNumThreads = 1;
static uint32_t const kDefaultWindowSize = 50;
static double const kDefaultPadeResolution = 10.0;
static double const kDefaultPadeMax = 1000.0;
static double const kDefaultMaxLkStart = 0.001;
static double const kDefaultMaxLkEnd = 0.3;
static double const kDefaultMaxLkResolution = 0.001;

CmdLineOptionsMaxLk::CmdLineOptionsMaxLk(std::string const base_command,
                                         int argc,
                                         char **argv,
                                         std::string version)
    : CmdLineOptions(version) {
  // Options allowed on command line and in configuration file.
  config_.add_options()
    ("num_threads",
     boost::program_options::value<uint32_t>
       (&num_threads_)->default_value(kDefaultNumThreads),
     "Number of threads.")
    ("seq_file,s",
     boost::program_options::value<std::string>(&seq_file_),
     "Sequence file.")
    ("lk_file,l",
     boost::program_options::value<std::string>
       (&lk_file_)->required(),
     "Two-site likelihood table.")
    ("pade_file,p",
     boost::program_options::value<std::string>(&pade_file_),
     "Pade coefficients.")
    ("prior_file,a",
     boost::program_options::value<std::string>(&prior_file_),
     "Prior on ancestral allele for each site.")
    ("mut_mat_file,m",
     boost::program_options::value<std::string>(&mut_mat_file_),
     "Mutation matrix.")
    ("window_size,w",
     boost::program_options::value<uint32_t>
       (&window_size_)->default_value(kDefaultWindowSize),
     "Window size.")
    ("pade_resolution",
     boost::program_options::value<double>
       (&pade_resolution_)->default_value(kDefaultPadeResolution),
     "Pade grid increment.")
    ("pade_max_rho",
     boost::program_options::value<double>(&pade_max_)
       ->default_value(kDefaultPadeMax),
     "Maximum Pade grid value.")
    ("pos_file",
     boost::program_options::value<std::string>(&pos_file_),
     "SNP positions for alternative input format.")
    ("snps_file",
     boost::program_options::value<std::string>(&snps_file_),
     "SNPs file for alternative input format.")
    ("max_lk_start",
     boost::program_options::value<double>(&max_lk_start_)
       ->default_value(kDefaultMaxLkStart),
     "Rho value to begin maximum likelihood estimation.")
    ("max_lk_end",
     boost::program_options::value<double>(&max_lk_end_)
       ->default_value(kDefaultMaxLkEnd),
     "Rho value to end maximum likelihood estimation.")
    ("max_lk_resolution",
     boost::program_options::value<double>(&max_lk_resolution_)
       ->default_value(kDefaultMaxLkResolution),
     "Amount to increment by for maximum likelihood estimation.");

  success_ = ParseOptions(base_command, argc, argv, version);
}

std::string CmdLineOptionsMaxLk::GetUsageString(
    std::string const base_command,
    std::string const program_name) const {
  return   "\nUsage: "
         + base_command
         + " "
         + program_name
         + " [options]\n";
}

bool CmdLineOptionsMaxLk::ProperInput(
  boost::program_options::variables_map const &variables_map) const {
  if (variables_map.count("seq_file") > 0 &&
      variables_map.count("pos_file") > 0 &&
      variables_map.count("snps_file") > 0) {
    fprintf(stderr,
            "Either a sequence file or a pair of SNP position "
            "and SNP sequence files should be used, but not "
            "both at the same time.\n");
    return false;
  }

  if (!(variables_map.count("seq_file") > 0 ||
        (variables_map.count("pos_file") > 0 &&
         variables_map.count("snps_file") > 0))) {
      fprintf(stderr, "Need sequence file.\n");
  }

  if (!(variables_map.count("lk_file") > 0)) {
      fprintf(stderr,
              "Need likeilhood table file.\n");
  }

  return variables_map["num_threads"].as<uint32_t>() >= 1 &&
         (variables_map.count("seq_file") > 0 ||
          (variables_map.count("pos_file") > 0 &&
           variables_map.count("snps_file") > 0)) &&
         variables_map.count("lk_file") > 0;
}
