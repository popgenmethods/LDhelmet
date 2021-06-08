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

#include "post_to_text/post_to_text_options.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <string>
#include <vector>

#include "common/command_line_options.h"

CmdLineOptionsPostToText::CmdLineOptionsPostToText(
    std::string const base_command,
    int argc,
    char **argv,
    std::string version) : CmdLineOptions(version) {
  // Options allowed on command line and in configuration file.
  config_.add_options()
    ("mean,m",
     boost::program_options::value(&mean_p_)
       ->zero_tokens()->default_value(false),
     "Specify option to output mean.")
    ("perc,p",
     boost::program_options::value<std::vector<double> >(&percs_)
       ->composing(),
     "Percentile value. Specify option multiple times for "
     "multiple percentiles."
     )
    ("output_file,o",
     boost::program_options::value<std::string>(&output_file_)
       ->required(),
     "Name of output file.");

  success_ = ParseOptions(base_command, argc, argv, version);
}

std::string CmdLineOptionsPostToText::GetUsageString(
    std::string const base_command,
    std::string const program_name) const {
  return "\nUsage: "
       + base_command
       + " "
       + program_name
       + " [options]"
       + " input_file\n";
}

bool CmdLineOptionsPostToText::ProperInput(
    boost::program_options::variables_map const &variables_map) const {
  if (variables_map.count("input_file") != 1) {
    fprintf(stderr, "Requires exactly one input file.\n");
  }

  return variables_map.count("input_file") == 1;
}
