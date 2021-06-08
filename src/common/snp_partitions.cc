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

#include "common/snp_partitions.h"

#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <cassert>
#include <string>
#include <utility>
#include <vector>

#include <boost/lexical_cast.hpp>

std::pair<SnpPartitions, SnpPartitions> ComputeSnpPartitions(
    std::vector<uint64_t> const &snp_pos,
    uint64_t const partition_length,
    uint64_t overlap) {
  if (partition_length < 2) {
    fprintf(stderr,
            "Partition length must be at least 2.\n");
    std::exit(1);
  }

  assert(snp_pos.size() >= 2);

  SnpPartitions snp_partitions;
  SnpPartitions snp_partitions_overlap;

  // Last partition takes up left over.
  assert(partition_length >= 2);
  uint64_t num_partitions = std::max(static_cast<uint64_t>(1),
                                     snp_pos.size()/(partition_length - 1));
  assert(num_partitions > 0);
  for (uint64_t partition_id = 0;
       partition_id < num_partitions - 1;
       ++partition_id) {
      assert(partition_id * partition_length < snp_pos.size());

      uint64_t start_snp_id = partition_id *(partition_length - 1);
      uint64_t end_snp_id = start_snp_id + partition_length;

      snp_partitions.push_back(std::make_pair(start_snp_id, end_snp_id));

      uint64_t start_snp_id_overlap = start_snp_id > overlap ?
                                        start_snp_id - overlap : 0;
      uint64_t end_snp_id_overlap =
        std::min(end_snp_id + overlap,
                 static_cast<uint64_t>(snp_pos.size()));
      snp_partitions_overlap.push_back(std::make_pair(start_snp_id_overlap,
                                                      end_snp_id_overlap));
  }

  // Last partititon.
  assert(num_partitions > 0);
  uint64_t last_start_snp_id = (num_partitions - 1) * (partition_length - 1);
  uint64_t last_end_snp_id = snp_pos.size();
  snp_partitions.push_back(std::make_pair(last_start_snp_id, last_end_snp_id));

  uint64_t last_start_snp_id_overlap =
    last_start_snp_id > overlap ?
      last_start_snp_id - overlap : 0;
  uint64_t last_end_snp_id_overlap = snp_pos.size();

  snp_partitions_overlap.push_back(std::make_pair(last_start_snp_id_overlap,
                                                  last_end_snp_id_overlap));

  assert(snp_partitions.size() == snp_partitions_overlap.size());
  assert(snp_partitions.size() > 0);
  assert(snp_partitions_overlap.size() > 0);

  return std::make_pair(snp_partitions, snp_partitions_overlap);
}

std::string ShowSnpPartitions(SnpPartitions const &snp_partitions) {
  std::string tmp;
  tmp += '[';
  for (int i = 0; i < static_cast<int>(snp_partitions.size()); ++i) {
    tmp += '[';
    tmp += boost::lexical_cast<std::string>(snp_partitions[i].first);
    tmp += ", ";
    tmp += boost::lexical_cast<std::string>(snp_partitions[i].second);
    tmp += ']';
  }
  tmp += ']';

  return tmp;
}
