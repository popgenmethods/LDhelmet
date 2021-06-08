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

#include "common/command_line_options.h"

#include <fstream>
#include <iostream>
#include <vector>

#include <boost/program_options.hpp>

std::string ShowCommandLineOptions(int argc, char **argv) {
  std::string tmp;

  assert(argc>0);
  tmp += std::string("Component: ") + argv[0] + '\n';
  tmp += std::string("Command line options:") + '\n';

  for (int i = 1; i < argc; ++i) {
      tmp += argv[i];
      tmp += ' ';
  }
  tmp += '\n';
  return tmp;
}

CmdLineOptions::CmdLineOptions(std::string const &version)
    : version_(version),
      success_(true),
      cmd_line_only_("General options"),
      config_("Component-specific options"),
      hidden_("Hidden options"),
      pos_options_() {
  // Options allowed only on the command line.
  cmd_line_only_.add_options()
    ("version,v", "Display version.")
    ("help,h", "Display help.");

  // Options allowed on command line or in configuration file but
  // are not shown to the user.
  hidden_.add_options()
    ("input_file",
     boost::program_options::value<std::vector<std::string> >
     (&input_files_),
     "Input file.")
    ("config,f",
     boost::program_options::value<std::string>(&config_file_),
     "Options file name.");

  // Positional arguments.
  pos_options_.add("input_file", -1);
}

bool CmdLineOptions::ParseOptions(std::string const base_command,
                                  int argc,
                                  char **argv,
                                  std::string const &version) {
  // Organize options into groups.
  boost::program_options::options_description cmd_line_options;
  cmd_line_options.add(cmd_line_only_).add(config_).add(hidden_);

  boost::program_options::options_description config_file_options;
  config_file_options.add(config_).add(hidden_);

  boost::program_options::options_description visible;
  visible.add(cmd_line_only_).add(config_);

  // Variables map (stores options).
  boost::program_options::variables_map variables_map;

  try {
    // Parse command line arguments.
    boost::program_options::store(
        boost::program_options::command_line_parser(argc, argv)
          .options(cmd_line_options).positional(pos_options_)
          .style(
              boost::program_options::command_line_style
                ::unix_style
              | boost::program_options::command_line_style
                ::allow_short)
          .run(),
        variables_map);

    // Parse options from configuration file.
    if (variables_map.count("config")) {
      std::string conf_file = variables_map["config"].as<std::string>();
      std::ifstream ifs(conf_file.c_str());
      if (!ifs) {
        fprintf(stderr, "Cannot open configuration file: %s.\n",
                config_file_.c_str());
        std::exit(1);
      } else {
        boost::program_options::store(
          boost::program_options::parse_config_file(ifs,
                                                    config_file_options),
          variables_map);
      }
    }

    // Special case handling for --version option.
    if (argc == 2 && variables_map.count("version")) {
      printf("%s\n", GetVersionString().c_str());
      return false;
    }

    // Special case handling for --help option.
    if (argc == 1 || variables_map.count("help")) {
      // Display version if requested.
      if (variables_map.count("version")) {
        printf("%s\n", GetVersionString().c_str());
      }

      // Display help text.
      printf("%s", GetUsageString(base_command, argv[0]).c_str());
      std::cout << visible;
      printf("\n");
      return false;
    }

    // Check for any exceptions from the command line parser.
    boost::program_options::notify(variables_map);

    // Bad options.
    if (!ProperInput(variables_map)) {
      // Display version if requested.
      if (variables_map.count("version")) {
        printf("%s\n", GetVersionString().c_str());
      }

      // Display help.
      printf("%s", GetUsageString(base_command, argv[0]).c_str());
      std::cout << visible;
      printf("\n");
      return false;
    }

    // Display version if user requests.
    if (variables_map.count("version")) {
      printf("%s\n", GetVersionString().c_str());
    }
  } catch (boost::program_options::required_option &e) {
    fprintf(stderr, "%s", GetUsageString(base_command, argv[0]).c_str());
    std::cerr << visible;
    fprintf(stderr, "\n");
    fprintf(stderr, "Command line parser: %s\n", e.what());
    return false;
  } catch (std::exception &e) {
    fprintf(stderr, "Command line parser: %s\n", e.what());
    return false;
  }

  return true;
}
