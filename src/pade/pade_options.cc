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

#include "pade/pade_options.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <string>

#include "common/command_line_options.h"

static uint32_t const kDefaultNumThreads = 1;
static uint64_t const kDefaultNumCoeffs = 11;
static double const kDefaultDefectThreshold = 40.0;

CmdLineOptionsPade::CmdLineOptionsPade(
    std::string const base_command,
    int argc,
    char **argv,
    std::string version)
    : CmdLineOptions(version) {
  // Options allowed on command line and in configuration file.
  config_.add_options()
    ("num_threads",
     boost::program_options::value<uint32_t>(&num_threads_)
       ->default_value(kDefaultNumThreads),
     "Number of threads to use.")
    ("conf_file,c",
     boost::program_options::value<std::string>(&conf_file_)
       ->required(),
     "Two-site configuration file.")
    ("output_file,o",
     boost::program_options::value<std::string>(&output_file_)
       ->required(),
     "Name for output file.")
    ("theta,t",
     boost::program_options::value<double>(&theta_)->required(),
     "Theta value.")
    ("num_coeff,x",
     boost::program_options::value<uint64_t>(&num_coeffs_)
       ->default_value(kDefaultNumCoeffs),
     "Number of Pade coefficients to compute.")
    ("defect_threshold",
     boost::program_options::value<double>(&defect_threshold_)
       ->default_value(kDefaultDefectThreshold),
     "Defect threshold for Pade coefficients.");

  success_ = ParseOptions(base_command, argc, argv, version);
}

std::string CmdLineOptionsPade::GetUsageString(
    std::string const base_command, std::string const program_name) const {
  return "\nUsage: "
       + base_command
       + " "
       + program_name
       + " [options]\n";
}

bool CmdLineOptionsPade::ProperInput(
  boost::program_options::variables_map const &variables_map) const {
  return variables_map["num_threads"].as<uint32_t>() >= 1 &&
         variables_map.count("conf_file") > 0 &&
         variables_map.count("output_file") > 0 &&
         variables_map.count("theta") > 0;
}
