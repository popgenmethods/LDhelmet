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

#ifndef LDHELMET_PADE_PADE_OPTIONS_H_
#define LDHELMET_PADE_PADE_OPTIONS_H_

// Header file implementing a basic command line and configuration file parser.

#include <string>

#include "common/command_line_options.h"

class CmdLineOptionsPade : public CmdLineOptions {
 public:
  CmdLineOptionsPade(std::string const base_command,
                     int argc,
                     char **argv,
                     std::string version);

  uint32_t num_threads_;

  uint64_t num_coeffs_;       // Number of pade coefficents to compute.
  std::string conf_file_;     // Configuration file for likelihood table.
  std::string output_file_;   // Name for output file.
  double theta_;              // Population-scaled mutation rate.
  double defect_threshold_;   // Defect threshold for pade coefficients.

 private:
  std::string GetUsageString(std::string const base_command,
                             std::string const program_name) const;

  bool ProperInput(
      boost::program_options::variables_map const &variables_map) const;
};

#endif  // LDHELMET_PADE_PADE_OPTIONS_H_
