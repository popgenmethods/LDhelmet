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

#include "common/read_confs.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include "common/conf.h"

// Load input configurations from file and construct conf index.
std::vector<Conf> LoadConfigurations(std::string const &conf_file_name) {
  FILE *fp = fopen(conf_file_name.c_str(), "r");
  if (fp == NULL) {
    fprintf(stderr,
            "Unable to open input configuration file: %s.\n",
            conf_file_name.c_str());
    std::exit(1);
  }
  std::vector<Conf> conf_list;

  size_t index = 0;
  while (!feof(fp)) {
    int tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    int number_read = fscanf(fp, "%d %d %d %d %d %d %d %d",
                             &tmp0, &tmp1, &tmp2, &tmp3,
                             &tmp4, &tmp5, &tmp6, &tmp7);
    if (number_read == 8) {
      assert(!feof(fp) && !ferror(fp));
      Conf conf;
      conf.a_ = tmp0;
      conf.A_ = tmp1;
      conf.b_ = tmp2;
      conf.B_ = tmp3;
      conf.ab_ = tmp4;
      conf.aB_ = tmp5;
      conf.Ab_ = tmp6;
      conf.AB_ = tmp7;
      assert(conf.ComputeDegree() >= 1);

      conf_list.push_back(conf);
      ++index;
    } else if (number_read != -1 || ferror(fp)) {
      fprintf(stderr, "Error parsing haplotype configuration file.\n");
      std::exit(1);
    } else {  // End of file.
      assert(feof(fp) && !ferror(fp));
    }
  }

  return conf_list;
}

boost::tuple<std::vector<size_t>, uint32_t, uint32_t, uint32_t>
PreprocessConfs(std::vector<Conf> const &conf_list) {
  // degree separation, max degree, max sample size, max locus
  boost::tuple<std::vector<size_t>, uint32_t, uint32_t, uint32_t>
    conf_meta_data;

  uint32_t old_degree = 0;
  conf_meta_data.get<0>().push_back(0);

  for (size_t conf_id = 0; conf_id < conf_list.size(); ++conf_id) {
    Conf const &conf = conf_list[conf_id];
    uint32_t degree = conf.ComputeDegree();
    if (degree > old_degree) {
      for (size_t i = old_degree; i < degree; ++i) {
        conf_meta_data.get<0>().push_back(conf_id);
      }
      old_degree = degree;
    }
    conf_meta_data.get<1>() = std::max(conf_meta_data.get<1>(), degree);
    conf_meta_data.get<2>() = std::max(conf_meta_data.get<2>(),
                                       conf.ComputeSize());
    conf_meta_data.get<3>() =
        std::max(conf_meta_data.get<3>(),
                 std::max(conf.ComputeaMar() + conf.ComputeAMar(),
                          conf.ComputebMar() + conf.ComputeBMar()));
  }
  conf_meta_data.get<0>().push_back(conf_list.size());

  assert(conf_meta_data.get<0>().size() == conf_meta_data.get<1>() + 2);
  return conf_meta_data;
}
