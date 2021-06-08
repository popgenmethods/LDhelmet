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

#ifndef LDHELMET_COMMON_LOAD_DATA_H_
#define LDHELMET_COMMON_LOAD_DATA_H_

#include <stdint.h>

#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include "common/mut_mat_prior.h"

boost::tuple<std::vector<std::string>, std::vector<uint64_t> > LoadSeqData(
    uint32_t const num_threads,
    std::string const &seq_file,
    std::string const &pos_file,
    std::string const &snps_file);

boost::tuple<MutationMatrix,
             std::vector<std::string>,
             std::vector<uint64_t>,
             Prior> LoadData(uint32_t const num_threads,
                             std::string const &mut_mat_file,
                             std::string const &seq_file,
                             std::string const &pos_file,
                             std::string const &snps_file,
                             std::string const &prior_file);

#endif  // LDHELMET_COMMON_LOAD_DATA_H_
