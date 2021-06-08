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

#include "find_confs/find_confs.h"

#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

#include "common/seq_file_parse.h"
#include "common/seq_process.h"
#include "common/version_number.h"
#include "find_confs/find_confs_options.h"

boost::tuple<std::vector<std::string>, std::vector<uint64_t> > FindSnpsHelper(
    std::vector<std::string> const &seqs,
    size_t begin,
    size_t end) {
  assert(seqs.size() > 0);
  size_t num_haps = seqs.size();

  std::vector<std::string> snp_seqs(num_haps);
  std::vector<uint64_t> snp_pos;

  for (size_t site = begin; site < end; ++site) {
    std::set<uint8_t> alleles;
    for (size_t hap = 0; hap < seqs.size(); ++hap) {
      uint8_t base = ConvertCharToUInt8(seqs[hap][site]);
      if (ValidSnpAlleleP(base)) {
        alleles.insert(base);
      }
    }

    if (alleles.size() == 2) {
      for (size_t hap = 0; hap < num_haps; ++hap) {
        snp_seqs[hap] += seqs[hap][site];
      }

      snp_pos.push_back(site);
    }
  }

  return boost::make_tuple(snp_seqs, snp_pos);
}

boost::tuple<std::vector<std::string>, std::vector<uint64_t> > FindSnps(
    uint32_t const num_threads,
    std::vector<std::string> const &seqs) {
  assert(seqs.size() > 0);

  std::vector<boost::tuple<std::vector<std::string>, std::vector<uint64_t> > >
    results(num_threads);

  // Partition sequence.
  size_t num_haps = seqs.size();
  size_t len_seq = seqs.front().size();
  size_t block_length = len_seq/num_threads;

  boost::thread_group tgroup;
  for (uint32_t thread = 0; thread < num_threads; ++thread) {
    size_t begin = thread * block_length;
    size_t end = thread < num_threads - 1 ?
                   (thread + 1) * block_length : len_seq;
    tgroup.create_thread(
        FindSnpsPartition(seqs, begin, end, results[thread]));
  }
  tgroup.join_all();

  uint64_t num_snps = results.front().get<0>().front().size();

  std::vector<std::string> snp_seqs(num_haps);
  std::vector<uint64_t> snp_pos;

  // Merge snp sequences.
  for (size_t hap = 0; hap < num_haps; ++hap) {
    snp_seqs.reserve(num_snps);
  }

  for (size_t hap = 0; hap < num_haps; ++hap) {
    for (size_t result_id = 0; result_id < results.size(); ++result_id) {
      snp_seqs[hap] += results[result_id].get<0>()[hap];
    }
  }

  // Merge snp positions.
  for (size_t result_id = 0; result_id < results.size(); ++result_id) {
    snp_pos.insert(snp_pos.end(),
                   results[result_id].get<1>().begin(),
                   results[result_id].get<1>().end());
  }

  assert(snp_seqs.front().size() == snp_pos.size());

  return boost::make_tuple(snp_seqs, snp_pos);
}

std::set<Conf> FindConfs(uint32_t const window_size,
                         std::vector<std::string> const &snp_seqs) {
  assert(snp_seqs.size() > 0);

  size_t len_seq = snp_seqs.front().size();
  std::set<Conf> confs;

  for (size_t site0 = 0; site0 < len_seq - 1; ++site0) {
    size_t max_site1 = std::min(site0 + window_size, len_seq);
    for (size_t site1 = site0 + 1; site1 < max_site1; ++site1) {
      Conf conf = DetermineConf(snp_seqs, site0, site1);

      // Original.
      confs.insert(conf);
      // Swap alleles of second locus.
      confs.insert(Conf(conf.a_,
                        conf.A_,
                        conf.B_,
                        conf.b_,
                        conf.aB_,
                        conf.ab_,
                        conf.AB_,
                        conf.Ab_));
      // Swap alleles of first locus.
      confs.insert(Conf(conf.A_,
                        conf.a_,
                        conf.b_,
                        conf.B_,
                        conf.Ab_,
                        conf.AB_,
                        conf.ab_,
                        conf.aB_));
      // Swap alleles at both loci.
      confs.insert(Conf(conf.A_,
                        conf.a_,
                        conf.B_,
                        conf.b_,
                        conf.AB_,
                        conf.Ab_,
                        conf.aB_,
                        conf.ab_));
    }
  }

  return confs;
}

std::set<Conf> MergeConfs(std::vector<std::set<Conf> > const &conf_sets) {
  std::set<Conf> confs;
  for (std::vector<std::set<Conf> >::const_iterator iter = conf_sets.begin();
    iter != conf_sets.end();
    ++iter) {
    confs.insert(iter->begin(), iter->end());
  }
  return confs;
}

std::set<Conf> FindConfsMT(uint32_t const num_threads,
                           uint32_t const window_size,
                           std::vector<std::string> const &snp_seqs) {
  size_t len_snps = snp_seqs.front().size();

  std::vector<std::set<Conf> > results(num_threads);

  // Partition sequences.
  size_t len_block = len_snps / num_threads;
  std::vector<size_t> breakpoints(num_threads + 1);

  for (size_t j = 0; j < breakpoints.size() - 1; ++j) {
    breakpoints[j] = j * len_block;
  }
  breakpoints.back() = len_snps;

  std::vector<std::vector<std::string> > new_snps(num_threads);
  for (size_t thread = 0; thread < new_snps.size(); ++thread) {
    new_snps[thread].resize(snp_seqs.size());
    size_t start = breakpoints[thread] > window_size + 1 ?
                     breakpoints[thread] - window_size + 1 : 0;
    size_t end = breakpoints[thread + 1];
    size_t len = end - start;
    for (size_t j = 0; j < snp_seqs.size(); ++j) {
      new_snps[thread][j] = snp_seqs[j].substr(start, len);
    }
  }

  boost::thread_group tgroup;
  for (uint32_t thread = 0; thread < num_threads; ++thread) {
    tgroup.create_thread(FindConfsPartition(window_size,
                                            new_snps[thread],
                                            results[thread]));
  }
  tgroup.join_all();

  return MergeConfs(results);
}

std::vector<Conf> FindConfsFromFiles(
    uint32_t const num_threads,
    uint32_t const window_size,
    std::vector<std::string> const &input_files) {
  std::vector<std::set<Conf> > final_results(input_files.size());

  // Check that every file can be opened before processing any of them.
  for (int i = 0; i < static_cast<int>(input_files.size()); ++i) {
    std::string const &input_file = input_files[i];
    FILE *fp = fopen(input_file.c_str(), "r");
    if (fp == NULL) {
      fprintf(stderr, "Could not open file %s.\n", input_file.c_str());
      std::exit(1);
    }
    fclose(fp);
  }

  for (size_t file_id = 0; file_id < input_files.size(); ++file_id) {
    std::vector<std::string> seqs = ReadSeqFile(input_files[file_id]);

    final_results[file_id] =
      FindConfsMT(num_threads,
                  window_size,
                  boost::get<0>(FindSnps(num_threads, seqs)));
  }

  std::set<Conf> confs_set = MergeConfs(final_results);
  std::vector<Conf> confs(confs_set.begin(), confs_set.end());
  return confs;
}
