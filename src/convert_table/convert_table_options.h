#ifndef LDHELMET_CONVERT_TABLE_CONVERT_TABLE_OPTIONS_H_
#define LDHELMET_CONVERT_TABLE_CONVERT_TABLE_OPTIONS_H_

#include <string>
#include <vector>

#include "common/command_line_options.h"

class CmdLineOptionsConvertTable : public CmdLineOptions {
public:
    CmdLineOptionsConvertTable(std::string const base_command,
                           int argc,
                           char **argv,
                           std::string version);
    
    uint32_t num_threads_;
    
    std::string input_file_;         // LDhat style .txt file containing all the likelihoods
    std::string output_file_;        // Name for output file.
    
    
private:
    std::string GetUsageString(std::string const base_command,
                               std::string const program_name) const;
    
    bool ProperInput(
                     boost::program_options::variables_map const &variables_map) const;
};

#endif  // LDHELMET_CONVERT_TABLE_CONVERT_TABLE_OPTIONS_H_
