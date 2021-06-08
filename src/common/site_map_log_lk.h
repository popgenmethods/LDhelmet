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

#ifndef LDHELMET_COMMON_SEQ_MAP_LOG_LK_H_
#define LDHELMET_COMMON_SEQ_MAP_LOG_LK_H_

#include <string>
#include <vector>

#include "common/rho_finder.h"
#include "common/lk_pade_table.h"
#include "common/mut_mat_prior.h"

class SiteMapLogLk {
 public:
  SiteMapLogLk(LkTable const &lk_table,
               uint32_t const window_size,
               std::vector<std::string> const &snp_seqs,
               uint64_t const start_snp_id,
               uint64_t const end_snmp_id,
               MutationMatrix const &mut_mat,
               Prior const &prior);

  uint32_t const window_size_;
  std::vector<double> site_log_lks_;
  RhoFinder rho_finder_;
  size_t len_snp_seq_;
  size_t rho_stride_;
};

#endif  // LDHELMET_COMMON_SEQ_MAP_LOG_LK_H_
