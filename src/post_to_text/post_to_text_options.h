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

#ifndef LDHELMET_POST_TO_TEXT_POST_TO_TEXT_OPTIONS_H_
#define LDHELMET_POST_TO_TEXT_POST_TO_TEXT_OPTIONS_H_

// Header file implementing a basic command line and configuration file parser.

#include <string>
#include <vector>

#include "common/command_line_options.h"

class CmdLineOptionsPostToText : public CmdLineOptions {
 public:
  CmdLineOptionsPostToText(std::string const base_command,
                           int argc,
                           char **argv,
                           std::string version);

  std::vector<double> percs_;  // Percentiles.
  bool mean_p_;                // Display mean.
  std::string output_file_;

 private:
  std::string GetUsageString(std::string const base_command,
                             std::string const program_name) const;

  bool ProperInput(
      boost::program_options::variables_map const &variables_map) const;
};

#endif  // LDHELMET_POST_TO_TEXT_POST_TO_TEXT_OPTIONS_H_
