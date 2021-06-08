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

#ifndef LDHELMET_COMMON_SEQ_FILE_PARSE_H_
#define LDHELMET_COMMON_SEQ_FILE_PARSE_H_

#include <stdint.h>

#include <string>
#include <vector>

std::vector<std::string> ReadFASTAFile(FILE *fp);

std::vector<std::string> ReadFASTQFile(FILE *fp);

std::vector<std::string> ReadSeqFile(std::string const &input_file);

std::vector<uint64_t> ReadPosFile(std::string const &pos_file);

#endif  // LDHELMET_COMMON_SEQ_FILE_PARSE_H_
