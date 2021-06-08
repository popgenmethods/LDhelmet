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

#ifndef LDHELMET_COMMON_SNP_PARTITIONS_H_
#define LDHELMET_COMMON_SNP_PARTITIONS_H_

#include <stdint.h>

#include <string>
#include <utility>
#include <vector>

typedef std::vector<std::pair<uint64_t, uint64_t> > SnpPartitions;

std::pair<SnpPartitions, SnpPartitions> ComputeSnpPartitions(
    std::vector<uint64_t> const &snp_pos,
    uint64_t const partition_length,
    uint64_t overlap);

std::string ShowSnpPartitions(SnpPartitions const &snp_partitions);

#endif  // LDHELMET_COMMON_SNP_PARTITIONS_H_
