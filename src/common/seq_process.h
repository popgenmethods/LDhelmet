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

#ifndef LDHELMET_COMMON_SEQ_PROCESS_H_
#define LDHELMET_COMMON_SEQ_PROCESS_H_

#include <cstdlib>
#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>

#include "common/conf.h"

uint8_t const kUndefAllele = 4;

typedef boost::array<uint8_t, 2> SnpAlleles;

inline bool ValidSnpAlleleP(uint8_t base) {
  return base < 4;
}

inline uint8_t ConvertCharToUInt8(char nucleotide) {
  uint8_t base;
  switch (nucleotide) {
    case 'A':
    case 'a':
      base = 0;
      break;
    case 'C':
    case 'c':
      base = 1;
      break;
    case 'G':
    case 'g':
      base = 2;
      break;
    case 'T':
    case 't':
      base = 3;
      break;
    default:
      base = 4;
      break;
  }
  return base;
}

// Does not assume site is diallelic.
SnpAlleles DetermineAllelesCareful(
  std::vector<std::string::const_iterator> &snp_iters);

// Assumes site is diallelic.
SnpAlleles DetermineAlleles(std::vector<std::string> const &snp_seqs,
                            size_t site);

Conf DetermineConf(std::vector<std::string> const &snp_seqs,
                   size_t site0,
                   size_t site1);

void CleanSnpAndPos(std::vector<std::string> *snp_seqs,
                    std::vector<uint64_t> *snp_pos);

#endif  // LDHELMET_COMMON_SEQ_PROCESS_H_
