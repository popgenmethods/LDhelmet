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

#include "common/site_map_log_lk.h"

#include <stdio.h>

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

#include <boost/array.hpp>

#include "common/lk_pade_table.h"
#include "common/mut_mat_prior.h"
#include "common/rho_finder.h"
#include "common/seq_process.h"

SiteMapLogLk::SiteMapLogLk(LkTable const &lk_table,
                           uint32_t const window_size,
                           std::vector<std::string> const &snp_seqs,
                           uint64_t const start_snp_id,
                           uint64_t const end_snp_id,
                           MutationMatrix const &mut_mat,
                           Prior const &prior)
    : window_size_(window_size), len_snp_seq_(), rho_stride_() {
  if (end_snp_id <= start_snp_id) {
    fprintf(stderr,
            "End SNP ID must be greater than start SNP ID.");
    std::exit(1);
  }

  assert(snp_seqs.size() >= 2);
  assert(snp_seqs.front().size() >= 2);

  ConfLkMap const &conf_lk_map = lk_table.conf_lk_map_;

  rho_finder_ = lk_table.rho_finder_;

  assert(snp_seqs.size() > 0);
  len_snp_seq_ = end_snp_id - start_snp_id;
  assert(len_snp_seq_ > 0);

  rho_stride_ = (len_snp_seq_ - 1) * window_size_;

  uint64_t size_site_log_lks = rho_stride_ * rho_finder_.rho_list_.size();

  if (size_site_log_lks > 2684354560) {
    fprintf(stderr,
            "Safety check: "
            "You probably don't want to analyze such a large number "
            "of SNPs at once. The amount of memory required will be "
            "at least 20 GB. If you do, you'll have to remove this "
            "safety check from the code.\n"
            "Recommendation: Use shorter partitions.\n");
    std::exit(1);
  }

  site_log_lks_ = std::vector<double>(size_site_log_lks);

  // Compute the likelihoods.
  for (size_t site0 = 0; site0 < len_snp_seq_ - 1; ++site0) {
    size_t max_site1 = std::min(site0 + window_size_, len_snp_seq_);
    for (size_t site1 = site0 + 1; site1 < max_site1; ++site1) {
      size_t offset_site0 = site0 + start_snp_id;
      size_t offset_site1 = site1 + start_snp_id;

      boost::array<Conf, 4> confs;
      // Original conf.
      confs[0] = DetermineConf(snp_seqs, offset_site0, offset_site1);
      // Swap alleles at second locus.
      confs[1] = Conf(confs[0].a_,
                      confs[0].A_,
                      confs[0].B_,
                      confs[0].b_,
                      confs[0].aB_,
                      confs[0].ab_,
                      confs[0].AB_,
                      confs[0].Ab_);
      // swap alleles at first locus
      confs[2] = Conf(confs[0].A_,
                      confs[0].a_,
                      confs[0].b_,
                      confs[0].B_,
                      confs[0].Ab_,
                      confs[0].AB_,
                      confs[0].ab_,
                      confs[0].aB_);
      // swap alleles at both loci
      confs[3] = Conf(confs[0].A_,
                      confs[0].a_,
                      confs[0].B_,
                      confs[0].b_,
                      confs[0].AB_,
                      confs[0].Ab_,
                      confs[0].aB_,
                      confs[0].ab_);

      SnpAlleles snp_alleles0 = DetermineAlleles(snp_seqs, offset_site0);
      SnpAlleles snp_alleles1 = DetermineAlleles(snp_seqs, offset_site1);

      for (size_t rho_id = 0;
           rho_id < rho_finder_.rho_list_.size();
           ++rho_id) {
        double lk =
            mut_mat[snp_alleles0[0]][snp_alleles0[1]]
              * mut_mat[snp_alleles1[0]][snp_alleles1[1]]
              * prior[offset_site0][snp_alleles0[0]]
              * prior[offset_site1][snp_alleles1[0]]
              * conf_lk_map.at(confs[0]).second[rho_id]
            + mut_mat[snp_alleles0[0]][snp_alleles0[1]]
              * mut_mat[snp_alleles1[1]][snp_alleles1[0]]
              * prior[offset_site0][snp_alleles0[0]]
              * prior[offset_site1][snp_alleles1[1]]
              * conf_lk_map.at(confs[1]).second[rho_id]
            + mut_mat[snp_alleles0[1]][snp_alleles0[0]]
              * mut_mat[snp_alleles1[0]][snp_alleles1[1]]
              * prior[offset_site0][snp_alleles0[1]]
              * prior[offset_site1][snp_alleles1[0]]
              * conf_lk_map.at(confs[2]).second[rho_id]
            + mut_mat[snp_alleles0[1]][snp_alleles0[0]]
              * mut_mat[snp_alleles1[1]][snp_alleles1[0]]
              * prior[offset_site0][snp_alleles0[1]]
              * prior[offset_site1][snp_alleles1[1]]
              * conf_lk_map.at(confs[3]).second[rho_id];

        site_log_lks_[rho_id * rho_stride_
                    + site0 * window_size_
                    + (site1 - site0)] = std::log(lk);
      }
    }
  }
}
