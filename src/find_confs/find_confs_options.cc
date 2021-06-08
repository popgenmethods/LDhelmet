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

#include "find_confs/find_confs_options.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <string>

static uint32_t const kDefaultNumThreads = 1;
static uint32_t const kDefaultWindowSize = 50;

CmdLineOptionsFindConfs::CmdLineOptionsFindConfs(std::string const base_command,
                                                 int argc,
                                                 char **argv,
                                                 std::string version)
    : CmdLineOptions(version) {
  // Options allowed on command line and in configuration file
  config_.add_options()
    ("num_threads",
     boost::program_options::value<uint32_t>(&num_threads_)
       ->default_value(kDefaultNumThreads),
     "Number of threads to use.")
    ("window_size,w",
     boost::program_options::value<uint32_t>(&window_size_)
       ->default_value(kDefaultWindowSize),
     "Window size.")
    ("output_file,o",
     boost::program_options::value<std::string>(&output_file_)
       ->required(),
     "Name for output file.");

  success_ = ParseOptions(base_command, argc, argv, version);
}


std::string CmdLineOptionsFindConfs::GetUsageString(
    std::string const base_command,
    std::string const program_name) const {
  return "\nUsage: "
       + base_command
       + " "
       + program_name
       + " [options] seq-file1 [seq-file2 ...]\n";
}

bool CmdLineOptionsFindConfs::ProperInput(
    boost::program_options::variables_map const &variables_map) const {
  if (variables_map.count("input_file") == 0) {
    fprintf(stderr, "Need at least one input file.\n");
  }

  if (variables_map.count("output_file") == 0) {
    fprintf(stderr, "Need output file name.\n");
  }

  return variables_map["num_threads"].as<uint32_t>() >= 1 &&
         variables_map.count("input_file")  > 0 &&
         variables_map.count("output_file") > 0;
}
