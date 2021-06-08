#include "convert_table/convert_table_options.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <string>
#include <vector>

#include "common/command_line_options.h"

CmdLineOptionsConvertTable::CmdLineOptionsConvertTable(std::string const base_command,
                                               int argc,
                                               char **argv,
                                               std::string version)
: CmdLineOptions(version) {
    config_.add_options()
    ("input_table",
     boost::program_options::value<std::string>(&input_file_)
     ->required(),
     "LDhat style table to be converted.")
    ("output_table",
     boost::program_options::value<std::string>(&output_file_)
     ->required(),
     "Name for output file.")
    ("config_file",
     boost::program_options::value<std::string>(&config_file_)
     ->default_value(""),
     "(Optional) File with configs.  This is only necessary if you have missing data.");
    
    success_ = ParseOptions(base_command, argc, argv, version);
}

std::string CmdLineOptionsConvertTable::GetUsageString(
                                                   std::string const base_command,
                                                   std::string const program_name) const {
    return "\nUsage: "
    + base_command
    + " "
    + program_name
    + " [options]\n";
}

bool CmdLineOptionsConvertTable::ProperInput(
                                         boost::program_options::variables_map const &variables_map) const {
    return variables_map.count("input_table") > 0 &&
    variables_map.count("output_table") > 0;
}
