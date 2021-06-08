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

#include "table_gen/output_writer.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <string>
#include <vector>

#include "common/conf.h"
#include "common/vector_definitions.h"
#include "table_gen/table_management.h"

InputConfBinaryWriter::InputConfBinaryWriter(
    std::string file_name,
    std::vector<Conf> const &conf_list,
    std::vector<size_t> const &degree_seps)
    : conf_list_(conf_list), degree_seps_(degree_seps) {
  fp_ = fopen(file_name.c_str(), "wb");
  if (fp_ == NULL) {
    fprintf(stderr,
            "Error opening file for InputConfBinaryWriter.\n");
    std::exit(1);
  }
}

void InputConfBinaryWriter::Write(uint32_t degree, Vec8 const &table) {
  for (size_t conf_id = degree_seps_[degree];
       conf_id < degree_seps_[degree + 1];
       ++conf_id) {
    // Write probability to file.
    double prob = GetTable(table, conf_list_[conf_id]);

    int num_written = fwrite(reinterpret_cast<const char *>(&prob),
                             sizeof(prob), 1, fp_);
    assert(num_written == 1);
  }
}
