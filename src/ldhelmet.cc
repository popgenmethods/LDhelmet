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

#include <stdio.h>

#include <cassert>
#include <cstdlib>
#include <utility>
#include <vector>

#include "common/version_number.h"
#include "find_confs/find_confs_component.h"
#include "max_lk/max_lk_component.h"
#include "pade/pade_component.h"
#include "post_to_text/post_to_text_component.h"
#include "rjmcmc/rjmcmc_component.h"
#include "table_gen/table_gen_component.h"
#include "convert_table/convert_table_component.h"

typedef int (*ComponentType)(std::string, int, char **);

int main(int argc, char **argv) {
  std::string version_string = VERSION_STRING;

  // List of components.
  std::vector<std::pair<std::string, ComponentType> >
      components;
  components.push_back(std::pair<std::string, ComponentType>
    ("find_confs", &FindConfsComponent));
  components.push_back(std::pair<std::string, ComponentType>
    ("table_gen", &TableGenComponent));
  components.push_back(std::pair<std::string, ComponentType>
    ("convert_table", &ConvertTableComponent));
  components.push_back(std::pair<std::string, ComponentType>
    ("pade", &PadeComponent));
  components.push_back(std::pair<std::string, ComponentType>
    ("rjmcmc", &RjmcmcComponent));
  components.push_back(std::pair<std::string, ComponentType>
    ("post_to_text", &PostToTextComponent));
  components.push_back(std::pair<std::string, ComponentType>
    ("max_lk", &MaxLkComponent));

  std::string command_list_str;
  for (std::vector<std::pair<std::string, ComponentType> >
         ::const_iterator citer = components.begin();
       citer != components.end();
       ++citer) {
    command_list_str += "  " + citer->first + "\n";
  }

  std::string help_str = "Version: "
                       + version_string
                       + "\n"
                       + std::string("Usage: ")
                       + argv[0]
                       + " <command> [args]\n\n"
                       + "Available commands:\n"
                       + command_list_str;

  if (argc < 2) {
    fprintf(stderr, "%s\n", help_str.c_str());
    std::exit(1);
  }

  std::string command = argv[1];

  if (command == "--help" || command == "-h") {
    printf("%s\n", help_str.c_str());
    std::exit(0);
  }

  if (command == "--version" || command == "-v") {
    printf("Version %s.\n", version_string.c_str());
    std::exit(0);
  }

  bool command_found = false;
  for (std::vector<std::pair<std::string, ComponentType> >
         ::const_iterator citer = components.begin();
       citer != components.end();
       ++citer) {
    if (command == citer->first) {
      command_found = true;
      std::string base_command = argv[0];
      (*citer->second)(base_command, argc - 1, argv + 1);
    }
  }

  if (!command_found) {
    fprintf(stderr, "\nCommand not found: %s\n\n", command.c_str());
    fprintf(stderr, "%s\n", help_str.c_str());
    std::exit(1);
  }
}
