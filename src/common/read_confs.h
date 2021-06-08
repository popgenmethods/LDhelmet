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

#ifndef LDHELMET_COMMON_READ_CONFS_INPUT_H_
#define LDHELMET_COMMON_READ_CONFS_INPUT_H_

#include <stddef.h>
#include <stdint.h>

#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include "common/conf.h"

std::vector<Conf> LoadConfigurations(std::string const &conf_file_name);

boost::tuple<std::vector<size_t>, uint32_t, uint32_t, uint32_t>
  PreprocessConfs(std::vector<Conf> const &conf_list);

std::vector<uint64_t> ReadPosFile(std::string const &pos_file);

#endif  // LDHELMET_COMMON_READ_CONFS_INPUT_H_
