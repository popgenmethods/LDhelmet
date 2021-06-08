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

#ifndef LDHELMET_COMMON_COMMAND_LINE_OPTIONS_H_
#define LDHELMET_COMMON_COMMAND_LINE_OPTIONS_H_

#include <string>
#include <vector>

#include <boost/program_options.hpp>

// Header file implementing a basic command line and configuration file parser.

// Returns user-specified command line arguments for display
std::string ShowCommandLineOptions(int argc, char **argv);

class CmdLineOptions {
 public:
  explicit CmdLineOptions(std::string const &version);

  virtual ~CmdLineOptions() { }

  inline bool success() const {
    return success_;
  }

  inline std::string GetVersionString() const {
    return "Version " + version_;
  }

  inline virtual std::string GetUsageString(
      std::string const base_command,
      std::string const program_name) const {
    return "<no help>";
  }

  virtual bool ProperInput(boost::program_options::variables_map
                             const &variables_map) const = 0;

  bool ParseOptions(std::string const base_command,
                    int argc,
                    char **argv,
                    std::string const &version);

  std::string version_;
  std::string config_file_;
  std::vector<std::string> input_files_;

 protected:
  bool success_;

  boost::program_options::options_description cmd_line_only_;
  boost::program_options::options_description config_;
  boost::program_options::options_description hidden_;
  boost::program_options::positional_options_description pos_options_;
};

#endif  // LDHELMET_COMMON_COMMAND_LINE_OPTIONS_H_
