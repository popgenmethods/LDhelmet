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

#ifndef LDHELMET_COMMON_HANDLE_OUTPUT_H_
#define LDHELMET_COMMON_HANDLE_OUTPUT_H_

#include <stdint.h>

#include <utility>
#include <vector>

#include "common/mut_mat_prior.h"
#include "common/snp_partitions.h"
#include "rjmcmc/rjmcmc_options.h"

std::vector<double> MergeResults(
    uint64_t num_snps,
    SnpPartitions const &snp_partitions,
    SnpPartitions const &snp_partitions_overlap,
    std::vector<std::vector<double> > const &results);

int64_t WriteHeaderToFile(FILE *fp,
                          uint64_t version_bit_string,
                          CmdLineOptionsRjmcmc const &cmd_line_options,
                          std::vector<uint64_t> const &snp_pos,
                          SnpPartitions const &snp_partitions,
                          SnpPartitions const &snp_partitions_overlap);

void WriteResultToFile(FILE *fp,
                       std::vector<uint64_t> const &snp_pos,
                       std::vector<double> const &result);

void WriteSamplesToFile(
    FILE *fp,
    std::vector<std::vector<std::pair<uint64_t, double> > >
      const &sample_store);

#endif  // LDHELMET_COMMON_HANDLE_OUTPUT_H_
