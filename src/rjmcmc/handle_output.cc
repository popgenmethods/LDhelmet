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

#include "rjmcmc/handle_output.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <utility>
#include <vector>
#include <algorithm>

#include <boost/functional/hash.hpp>

#include "common/mut_mat_prior.h"
#include "common/snp_partitions.h"
#include "rjmcmc/rjmcmc_options.h"

std::vector<double> MergeResults(
    uint64_t num_snps,
    SnpPartitions const &snp_partitions,
    SnpPartitions const &snp_partitions_overlap,
    std::vector<std::vector<double> > const &results) {
  std::vector<double> rho_map(num_snps - 1);
  for (size_t partition_id = 0;
       partition_id < snp_partitions_overlap.size();
       ++partition_id) {
    assert(snp_partitions[partition_id].first >=
             snp_partitions_overlap[partition_id].first);
    assert(snp_partitions[partition_id].second >
             snp_partitions_overlap[partition_id].first);

    size_t start_copy = snp_partitions[partition_id].first
                       - snp_partitions_overlap[partition_id].first;
    size_t end_copy = snp_partitions[partition_id].second
                     - snp_partitions_overlap[partition_id].first;

    assert(end_copy > start_copy);
    assert(start_copy < results[partition_id].size());
    assert(end_copy > 0);
    assert(end_copy - 1 <= results[partition_id].size());

    assert(!((end_copy - 1) - start_copy + snp_partitions[partition_id].first
             > rho_map.size()));

    std::copy(results[partition_id].begin() + start_copy,
              results[partition_id].begin() + end_copy - 1,
              rho_map.begin() + snp_partitions[partition_id].first);
  }
  return rho_map;
}

int64_t WriteHeaderToFile(FILE *fp,
                                 uint64_t version_bit_string,
                                 CmdLineOptionsRjmcmc const &cmd_line_options,
                                 std::vector<uint64_t> const &snp_pos,
                                 SnpPartitions const &snp_partitions,
                                 SnpPartitions const &snp_partitions_overlap) {
  int num_written;
  num_written = fwrite(reinterpret_cast<char const *>(&version_bit_string),
                       sizeof(version_bit_string), 1, fp);
  assert(num_written == 1);

  uint64_t seed_to_file = static_cast<uint64_t>(cmd_line_options.seed_);
  num_written = fwrite(reinterpret_cast<char const *>(&seed_to_file),
                       sizeof(seed_to_file), 1, fp);
  assert(num_written == 1);

  uint64_t burn_in_to_file = cmd_line_options.burn_in_;
  num_written = fwrite(reinterpret_cast<char const *>(&burn_in_to_file),
                       sizeof(burn_in_to_file), 1, fp);
  assert(num_written == 1);

  uint64_t num_iter_to_file = cmd_line_options.num_iter_;
  num_written = fwrite(reinterpret_cast<char const *>(&num_iter_to_file),
                       sizeof(num_iter_to_file), 1, fp);
  assert(num_written == 1);

  uint64_t stats_thin_to_file = cmd_line_options.stats_thin_;
  num_written = fwrite(reinterpret_cast<char const *>(&stats_thin_to_file),
                       sizeof(stats_thin_to_file), 1, fp);
  assert(num_written == 1);

  double block_penalty_to_file = cmd_line_options.block_penalty_;
  num_written = fwrite(reinterpret_cast<char const *>(&block_penalty_to_file),
                       sizeof(block_penalty_to_file), 1, fp);
  assert(num_written == 1);

  double prior_rate_mean_to_file = cmd_line_options.prior_rate_mean_;
  num_written = fwrite(reinterpret_cast<char const *>(&prior_rate_mean_to_file),
                       sizeof(prior_rate_mean_to_file), 1, fp);
  assert(num_written == 1);

  uint64_t window_size_to_file = cmd_line_options.window_size_;
  num_written = fwrite(reinterpret_cast<char const *>(&window_size_to_file),
                       sizeof(window_size_to_file), 1, fp);
  assert(num_written == 1);

  uint64_t partition_length_to_file = cmd_line_options.partition_length_;
  num_written = fwrite(reinterpret_cast<char const *>(
                           &partition_length_to_file),
                       sizeof(partition_length_to_file), 1, fp);
  assert(num_written == 1);

  uint64_t overlap_length_to_file = cmd_line_options.overlap_length_;
  num_written = fwrite(reinterpret_cast<char const *>(&overlap_length_to_file),
                       sizeof(overlap_length_to_file), 1, fp);
  assert(num_written == 1);

  double pade_resolution_to_file = cmd_line_options.pade_resolution_;
  num_written = fwrite(reinterpret_cast<char const *>(&pade_resolution_to_file),
                       sizeof(pade_resolution_to_file), 1, fp);
  assert(num_written == 1);

  double pade_max_to_file = cmd_line_options.pade_max_;
  num_written = fwrite(reinterpret_cast<char const *>(&pade_max_to_file),
                       sizeof(pade_max_to_file), 1, fp);
  assert(num_written == 1);

  double max_lk_start_to_file = cmd_line_options.max_lk_start_;
  num_written = fwrite(reinterpret_cast<char const *>(&max_lk_start_to_file),
                       sizeof(max_lk_start_to_file), 1, fp);
  assert(num_written == 1);

  double max_lk_end_to_file = cmd_line_options.max_lk_end_;
  num_written = fwrite(reinterpret_cast<char const *>(&max_lk_end_to_file),
                       sizeof(max_lk_end_to_file), 1, fp);
  assert(num_written == 1);

  double max_lk_resolution_to_file = cmd_line_options.max_lk_resolution_;
  num_written = fwrite(reinterpret_cast<char const *>(
                           &max_lk_resolution_to_file),
                       sizeof(max_lk_resolution_to_file), 1, fp);
  assert(num_written == 1);

  uint64_t num_snps_to_file = snp_pos.size();
  num_written = fwrite(reinterpret_cast<char const *>(&num_snps_to_file),
                       sizeof(num_snps_to_file), 1, fp);
  assert(num_written == 1);

  // Write snp positions to file.
  for (size_t i = 0; i < snp_pos.size(); ++i) {
    uint64_t snp_pos1 = snp_pos[i];
    num_written = fwrite(reinterpret_cast<char const *>(&snp_pos1),
                         sizeof(snp_pos1), 1, fp);
    assert(num_written == 1);
  }

  uint64_t num_snp_partitions_to_file = snp_partitions.size();
  num_written = fwrite(reinterpret_cast<char const *>(
                           &num_snp_partitions_to_file),
                       sizeof(num_snp_partitions_to_file), 1, fp);
  assert(num_written == 1);

  // Write partition ranges to file.
  for (size_t i = 0; i < snp_partitions.size(); ++i) {
    uint64_t begin = snp_partitions[i].first;
    uint64_t end = snp_partitions[i].second;
    num_written = fwrite(reinterpret_cast<char const *>(&begin),
                         sizeof(begin), 1, fp);
    assert(num_written == 1);
    num_written = fwrite(reinterpret_cast<char const *>(&end),
                         sizeof(end), 1, fp);
    assert(num_written == 1);
  }

  assert(snp_partitions.size() == snp_partitions_overlap.size());
  // Write partition ranges with overlap to file.
  for (size_t i = 0; i < snp_partitions_overlap.size(); ++i) {
    uint64_t begin = snp_partitions_overlap[i].first;
    uint64_t end = snp_partitions_overlap[i].second;
    num_written = fwrite(reinterpret_cast<char const *>(&begin),
                         sizeof(begin), 1, fp);
    assert(num_written == 1);
    num_written = fwrite(reinterpret_cast<char const *>(&end),
                         sizeof(end), 1, fp);
    assert(num_written == 1);
  }

  // Create room in the file for the marginal expectations.
  int64_t marg_exp_pos = ftell(fp);
  if (marg_exp_pos == -1) {
    fprintf(stderr, "Error writing to output file.\n");
    std::exit(1);
  }

  double nil = 0.0;
  for (size_t i = 0; i < snp_pos.size(); ++i) {
    num_written = fwrite(reinterpret_cast<char const *>(&nil),
                         sizeof(nil), 1, fp);
    assert(num_written == 1);
  }
  fflush(fp);

  return marg_exp_pos;
}

void WriteResultToFile(FILE *fp,
                       std::vector<uint64_t> const &snp_pos,
                       std::vector<double> const &result) {
  if (fp == NULL) {
    fprintf(stderr, "Error opening output file to write result.\n");
    std::exit(1);
  }

  assert(snp_pos.size() == result.size() + 1);

  for (size_t i = 0; i < result.size(); ++i) {
    double rate = result.at(i);
    int num_written = fwrite(reinterpret_cast<char *>(&rate),
                             sizeof(rate), 1, fp);
    assert(num_written == 1);
  }

  {
    double rate = -1.0;
    int num_written = fwrite(reinterpret_cast<char const *>(&rate),
                             sizeof(rate), 1, fp);
    assert(num_written == 1);
  }
  fflush(fp);
}

void WriteSamplesToFile(
    FILE *fp,
    std::vector<std::vector<std::pair<uint64_t, double> > >
      const &sample_store) {
  size_t seed = 0;

  uint64_t num_samps = sample_store.size();
  boost::hash_combine(seed, num_samps);

  for (size_t i = 0; i < sample_store.size(); ++i) {
    for (size_t j = 0; j < sample_store[i].size(); ++j) {
      uint64_t snp_pos1 = sample_store[i][j].first;
      double rate = sample_store[i][j].second;
      boost::hash_combine(seed, snp_pos1);
      boost::hash_combine(seed, rate);
    }
  }

  uint64_t seed_to_file = static_cast<uint64_t>(seed);
  int num_written;
  num_written = fwrite(reinterpret_cast<char const *>(&seed_to_file),
                       sizeof(seed_to_file), 1, fp);
  assert(num_written == 1);
  num_written = fwrite(reinterpret_cast<char const *>(&num_samps),
                       sizeof(num_samps), 1, fp);
  assert(num_written == 1);

  for (size_t i = 0; i < sample_store.size(); ++i) {
    uint64_t num_change_points = sample_store[i].size();

    num_written = fwrite(reinterpret_cast<char const *>(&num_change_points),
                         sizeof(num_change_points), 1, fp);
    for (size_t j = 0; j < sample_store[i].size(); ++j) {
      uint64_t snp_pos1 = sample_store[i][j].first;
      double rate = sample_store[i][j].second;
      num_written = fwrite(reinterpret_cast<char const *>(&snp_pos1),
                           sizeof(snp_pos1), 1, fp);
      assert(num_written == 1);
      num_written = fwrite(reinterpret_cast<char const *>(&rate),
                           sizeof(rate), 1, fp);
      assert(num_written == 1);
    }
  }
}
