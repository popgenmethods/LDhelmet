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

#ifndef LDHELMET_COMMON_MUT_MAT_PRIOR_H_
#define LDHELMET_COMMON_MUT_MAT_PRIOR_H_

#include <stdint.h>

#include <string>
#include <vector>

#include <boost/array.hpp>

typedef boost::array<boost::array<double, 4>, 4> MutationMatrix;

typedef boost::array<double, 4> PriorSite;
typedef std::vector<PriorSite> Prior;

MutationMatrix NormalizedMutMatrix(MutationMatrix const &mut_mat);

MutationMatrix LoadMutationMatrix(std::string const &mut_mat_file);

MutationMatrix DefaultMutMatrix();

Prior LoadPrior(MutationMatrix const &mut_mat,
                std::vector<uint64_t> const &snp_pos,
                std::string const &prior_file);

Prior DefaultPrior(MutationMatrix const &mut_mat,
                   std::vector<uint64_t> const &snp_pos);

#endif  // LDHELMET_COMMON_MUT_MAT_PRIOR_H_
