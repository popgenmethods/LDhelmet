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

#include "find_confs/find_confs_component.h"

#include <stdint.h>
#include <stdio.h>

#include <vector>

#include "common/conf.h"
#include "common/version_number.h"
#include "find_confs/find_confs.h"
#include "find_confs/find_confs_options.h"

int FindConfsComponent(std::string const base_command,
                       int argc,
                       char **argv) {
  std::string version_string = FIND_CONFS_VERSION_STRING;

  CmdLineOptionsFindConfs cmd_line_options(base_command,
                                           argc,
                                           argv,
                                           version_string);

  if (!cmd_line_options.success()) {
    std::exit(1);
  }

  printf("%s\n", ShowCommandLineOptions(argc, argv).c_str());

  uint32_t const num_threads = cmd_line_options.num_threads_;
  uint32_t const window_size = cmd_line_options.window_size_;

  printf("Number of threads: %d.\n", static_cast<int>(num_threads));
  assert(cmd_line_options.num_threads_ > 0);
  printf("Window size: %d.\n", static_cast<int>(window_size));

  if (window_size < 2) {
    fprintf(stderr, "Window size must be at least 2.\n");
    std::exit(1);
  }

  std::string const &output_file = cmd_line_options.output_file_;
  std::vector<std::string> const &input_files = cmd_line_options.input_files_;

  std::vector<Conf> confs =
    FindConfsFromFiles(num_threads,
                       window_size,
                       input_files);

  printf("Number of configurations: %d.\n", static_cast<int>(confs.size()));

  {
    printf("Writing configurations to %s.\n", output_file.c_str());

    FILE *fp = fopen(output_file.c_str(), "w");
    if (fp == NULL) {
      fprintf(stderr, "Error opening output file for writing "
                      "haplotype configurations.\n");
      std::exit(1);
    }
    for (int i = 0; i < static_cast<int>(confs.size()); ++i) {
      fprintf(fp, "%s\n", confs[i].GetString().c_str());
    }
    assert(fp != NULL);
    fclose(fp);
  }

  return 0;
}
