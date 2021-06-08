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

#ifndef LDHELMET_COMMON_LK_PADE_TABLE_H_
#define LDHELMET_COMMON_LK_PADE_TABLE_H_

#include <stdint.h>

#include <string>
#include <utility>
#include <vector>

#include <boost/unordered_map.hpp>

#include "common/conf.h"
#include "common/rho_finder.h"

class PadeCoefficients {
 public:
  PadeCoefficients() { }
  explicit PadeCoefficients(size_t num_coeffs) : coeffs_(num_coeffs) { }

  std::vector<double> coeffs_;
  std::vector<std::vector<double> > roots_;
};

typedef boost::unordered_map<
  Conf, std::pair<uint64_t, std::vector<double> >, ConfHash> ConfLkMap;
typedef boost::unordered_map<
  Conf, std::pair<uint64_t, PadeCoefficients>, ConfHash> ConfPadeMap;

class LkTable {
 public:
  LkTable();

  uint64_t version_bit_string_;
  uint64_t version_number_;
  uint64_t num_confs_;
  double theta_;
  uint8_t interpolate_;

  uint64_t num_lk_rhos_;
  ConfLkMap conf_lk_map_;
  RhoFinder rho_finder_;
};

class PadeTable {
 public:
  PadeTable();

  uint64_t version_bit_string_;
  uint64_t version_number_;
  uint64_t num_confs_;
  double theta_;
  uint64_t num_coeffs_;
  double defect_threshold_;

  ConfPadeMap conf_pade_map_;
};

class PadeComputation {
 public:
  PadeComputation(LkTable *lk_table,
                  PadeTable const &pade_table,
                  ConfPadeMap::const_iterator start_iter,
                  ConfPadeMap::const_iterator end_iter);

  double ComputePadeLk(PadeCoefficients const &pade_coefficients,
                       double pade_defect_threshold,
                       double rho);

  void operator()();

  LkTable *lk_table_;
  PadeTable const &pade_table_;
  ConfLkMap &conf_lk_map_;
  ConfPadeMap const &conf_pade_map_;
  ConfPadeMap::const_iterator start_iter_, end_iter_;
};

LkTable LoadLikelihoodAndPade(std::string const &lk_file,
                              std::string const &pade_file,
                              double const pade_resolution,
                              double const pade_max,
                              uint32_t const window_size,
                              std::vector<std::string> const &snp_seqs,
                              uint32_t const num_threads);

#endif  // LDHELMET_COMMON_LK_PADE_TABLE_H_
