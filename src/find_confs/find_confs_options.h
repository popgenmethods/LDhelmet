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

#ifndef LDHELMET_FIND_CONFS_FIND_CONFS_OPTIONS
#define LDHELMET_FIND_CONFS_FIND_CONFS_OPTIONS

#include <string>

#include "common/command_line_options.h"

class CmdLineOptionsFindConfs : public CmdLineOptions {
 public:
  CmdLineOptionsFindConfs(std::string const base_command,
                          int argc,
                          char **argv,
                          std::string version);

  uint32_t num_threads_;
  uint32_t window_size_;

  std::string output_file_;  // Output sample configurations file.

 private:
  std::string GetUsageString(std::string const base_command,
                             std::string const program_name) const;

  // Predicate to determine if options are proper.
  bool ProperInput(
      boost::program_options::variables_map const &variables_map) const;
};

#endif  // LDHELMET_FIND_CONFS_FIND_CONFS_OPTIONS
