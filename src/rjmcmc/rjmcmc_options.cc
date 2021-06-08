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

#include "rjmcmc/rjmcmc_options.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <string>

#include "common/command_line_options.h"
#include "rjmcmc/ran_num_gen.h"

static uint32_t const kDefaultNumThreads = 1;
static uint64_t const kDefaultStatsThin = 1000;
static uint32_t const kDefaultNumIter = 10000;
static uint32_t const kDefaultBurnIn = 1000;
static uint32_t const kDefaultWindowSize = 50;
static double const kDefaultBlockPenalty = 10.0;
static double const kDefaultPriorRateMean = 1.0;
static double const kDefaultPadeResolution = 10;
static double const kDefaultPadeMax = 1000;
static uint64_t const kDefaultPartitionLength = 4001;
static uint64_t const kDefaultOverlapLength = 200;
static double const kDefaultMaxLkStart = 0.001;
static double const kDefaultMaxLkEnd = 0.3;
static double const kDefaultMaxLkResolution = 0.001;

CmdLineOptionsRjmcmc::CmdLineOptionsRjmcmc(std::string const base_command,
                                           int argc,
                                           char **argv,
                                           std::string const &version)
    : CmdLineOptions(version) {
  // Options allowed on command line and in configuration file.
  config_.add_options()
    ("seed",
     boost::program_options::value<RanNumGen::SeedType>(&seed_)
       ->default_value(kDefaultSeed),
     "Seed for pseudo-random number generator.")
    ("num_threads",
     boost::program_options::value<uint32_t>(&num_threads_)
       ->default_value(kDefaultNumThreads),
     "Number of threads.")
    ("output_file,o",
     boost::program_options::value<std::string>(&output_file_)
       ->required(),
     "Name of output file.")
    ("stats_thin",
     boost::program_options::value<uint64_t>(&stats_thin_)
       ->default_value(kDefaultStatsThin),
     "Thinning parameter for summary statistics.")
    ("num_iter,n",
     boost::program_options::value<uint32_t>(&num_iter_)
       ->default_value(kDefaultNumIter),
     "Number of iterations to run rjMCMC.")
    ("burn_in",
     boost::program_options::value<uint32_t>(&burn_in_)
       ->default_value(kDefaultBurnIn),
     "Number of iterations for burn-in "
     "(in addition to number of iterations to run rjMCMC).")
    ("block_penalty,b",
     boost::program_options::value<double>(&block_penalty_)
       ->default_value(kDefaultBlockPenalty),
     "Block penalty for rjMCMC.")
    ("prior_rate",
     boost::program_options::value<double>(&prior_rate_mean_)
       ->default_value(kDefaultPriorRateMean),
     "Prior mean on recombination rate.")
    ("seq_file,s",
     boost::program_options::value<std::string>(&seq_file_),
     "Sequence file.")
    ("lk_file,l",
     boost::program_options::value<std::string>(&lk_file_)
       ->required(),
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
     boost::program_options::value<uint32_t>(&window_size_)
       ->default_value(kDefaultWindowSize),
     "Window size.")
    ("pade_resolution",
     boost::program_options::value<double>(&pade_resolution_)
       ->default_value(kDefaultPadeResolution),
     "Pade grid increment.")
    ("pade_max_rho",
     boost::program_options::value<double>(&pade_max_)
       ->default_value(kDefaultPadeMax),
     "Maximum Pade grid value.")
    ("partition_length",
     boost::program_options::value<uint64_t>(&partition_length_)
       ->default_value(kDefaultPartitionLength),
     "Partition length (number of SNPs).")
    ("overlap_length",
     boost::program_options::value<uint64_t>(&overlap_length_)
       ->default_value(kDefaultOverlapLength),
     "Overlap length.")
    ("pos_file",
     boost::program_options::value<std::string>(&pos_file_),
     "SNP positions for alternative input format.")
    ("snps_file",
     boost::program_options::value<std::string>(&snps_file_),
     "SNPs file for alternative input format.")
    ("max_lk_start",
     boost::program_options::value<double>(&max_lk_start_)
       ->default_value(kDefaultMaxLkStart),
     "Rho value to begin maximum likelihood estimation of "
     "background rate.")
    ("max_lk_end",
     boost::program_options::value<double>(&max_lk_end_)
       ->default_value(kDefaultMaxLkEnd),
     "Rho value to end maximum likelihood estimation of "
     "background rate.")
    ("max_lk_resolution",
     boost::program_options::value<double>(&max_lk_resolution_)
       ->default_value(kDefaultMaxLkResolution),
     "Amount to increment by for maximum likelihood estimation "
     "of background rate.");

  success_ = ParseOptions(base_command, argc, argv, version);
}

std::string CmdLineOptionsRjmcmc::GetUsageString(
    std::string const base_command, std::string const program_name) const {
  return "\nUsage: "
       + base_command
       + " "
       + program_name
       + " [options]\n";
}

bool CmdLineOptionsRjmcmc::ProperInput(
    boost::program_options::variables_map const &variablesMap) const {
  if (variablesMap.count("seq_file") > 0 &&
      variablesMap.count("pos_file") > 0 &&
      variablesMap.count("snps_file") > 0) {
    fprintf(stderr,
            "Either a sequence file or a pair of SNP position "
            "and SNP sequence files should be used, but not "
            "both at the same time.\n");
    return false;
  }

  return variablesMap["num_threads"].as<uint32_t>() >= 1 &&
         variablesMap.count("output_file") > 0 &&
         (variablesMap.count("seq_file") > 0 ||
          (variablesMap.count("pos_file") > 0 &&
              variablesMap.count("snps_file") > 0)) &&
         variablesMap.count("lk_file") > 0;
}
