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

#include "common/load_data.h"

#include <stdio.h>

#include <cassert>
#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include "common/mut_mat_prior.h"
#include "common/seq_file_parse.h"
#include "common/seq_process.h"
#include "find_confs/find_confs.h"

boost::tuple<std::vector<std::string>, std::vector<uint64_t> > LoadSeqData(
    uint32_t const num_threads,
    std::string const &seq_file,
    std::string const &pos_file,
    std::string const &snps_file) {
  std::vector<std::string> snp_seqs;
  std::vector<uint64_t> snp_pos;

  if (!seq_file.empty()) {
    printf("Finding SNPs in sequence file.\n");
    boost::tie(snp_seqs, snp_pos) = FindSnps(num_threads,
                                             ReadSeqFile(seq_file));
  } else if (!pos_file.empty() &&
            !snps_file.empty()) {
    printf("Loading SNPs.\n");

    snp_seqs = ReadSeqFile(snps_file);

    printf("Loading SNP positions.\n");

    snp_pos = ReadPosFile(pos_file);

    if (snp_seqs.front().size() != snp_pos.size()) {
      fprintf(stderr,
              "Number of SNPs in SNP sequence file is different from number "
              "of positions in SNP position file.\n");
      std::exit(1);
    }

    // Remove non-diallelic snps from both snp_seqs and snp_pos.
    // Note: Modifies in place.
    CleanSnpAndPos(&snp_seqs, &snp_pos);
  } else {
    fprintf(stderr,
           "No input files found. The command-line parser should have "
           "detected this earlier, so there's a bug in the code.\n");
    std::exit(1);
  }

  printf("Number of haplotypes: %d.\n", static_cast<int>(snp_seqs.size()));
  printf("Number of SNPs: %d.\n", static_cast<int>(snp_pos.size()));

  assert(snp_seqs.front().size() == snp_pos.size());

  if (snp_seqs.size() < 2 || snp_pos.size() < 2) {
    fprintf(stderr,
            "Number of sequences must be at least 2 and number of "
            "SNP positions must be at least 2.\n");
    std::exit(1);
  }

  return boost::make_tuple(snp_seqs, snp_pos);
}

boost::tuple<MutationMatrix,
             std::vector<std::string>,
             std::vector<uint64_t>,
             Prior> LoadData(uint32_t const num_threads,
                             std::string const &mut_mat_file,
                             std::string const &seq_file,
                             std::string const &pos_file,
                             std::string const &snps_file,
                             std::string const &prior_file) {
  MutationMatrix mut_mat = LoadMutationMatrix(mut_mat_file);

  std::vector<std::string> snp_seqs;
  std::vector<uint64_t> snp_pos;

  boost::tie(snp_seqs, snp_pos) = LoadSeqData(num_threads,
                                              seq_file,
                                              pos_file,
                                              snps_file);

  Prior prior = LoadPrior(mut_mat, snp_pos, prior_file);

  return boost::make_tuple(mut_mat, snp_seqs, snp_pos, prior);
}
