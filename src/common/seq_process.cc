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

#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <string>
#include <vector>

#include "common/seq_process.h"

SnpAlleles DetermineAllelesCareful(
  std::vector<std::string::const_iterator> &snp_iters) {
  SnpAlleles snp_alleles = {{kUndefAllele, kUndefAllele}};  // {'N','N'}

  size_t hap = 0;

  // Find first allele.
  for (; hap < snp_iters.size(); ++hap) {
    uint8_t base = ConvertCharToUInt8(*snp_iters[hap]);
    if (ValidSnpAlleleP(base)) {
      snp_alleles[0] = base;
      ++hap;
      break;
    }
  }
  // Find second allele.
  for (; hap < snp_iters.size(); ++hap) {
    uint8_t base = ConvertCharToUInt8(*snp_iters[hap]);
    if (ValidSnpAlleleP(base) && base != snp_alleles[0]) {
      snp_alleles[1] = base;
      ++hap;
      break;
    }
  }

  // Check if site is diallelic.
  for (; hap < snp_iters.size(); ++hap) {
    uint8_t base = ConvertCharToUInt8(*snp_iters[hap]);
    if (!(base == kUndefAllele  ||
       base == snp_alleles[0] ||
       base == snp_alleles[1])) {
      // Not diallelic; set snp_alleles to {N,N} to indicate an
      // invalid site.
      snp_alleles[0] = kUndefAllele;
      snp_alleles[1] = kUndefAllele;
      break;
    }
  }

  // Order alleles according to (A,C,G,N,T) (alphabetical).
  if (snp_alleles[0] > snp_alleles[1]) {
      std::swap(snp_alleles[0], snp_alleles[1]);
  }

  return snp_alleles;
}

SnpAlleles DetermineAlleles(std::vector<std::string> const &snp_seqs,
                            size_t site) {
  SnpAlleles snp_alleles = {{kUndefAllele, kUndefAllele}};  // {'N','N'}

  size_t hap = 0;

  // Find first allele.
  for (; hap < snp_seqs.size(); ++hap) {
    uint8_t base = ConvertCharToUInt8(snp_seqs[hap][site]);
    if (ValidSnpAlleleP(base)) {
      snp_alleles[0] = base;
      ++hap;
      break;
    }
  }
  // Find second allele.
  for (; hap < snp_seqs.size(); ++hap) {
    uint8_t base = ConvertCharToUInt8(snp_seqs[hap][site]);
    if (ValidSnpAlleleP(base) && base != snp_alleles[0]) {
      snp_alleles[1] = base;
      ++hap;
      break;
    }
  }

  if (!ValidSnpAlleleP(snp_alleles[0]) || !ValidSnpAlleleP(snp_alleles[1])) {
    fprintf(stderr,
            "A SNP site has been detected as non-diallelic. "
            "This should have been caught earlier. "
            "There is probably a bug in the code.\n");
    std::exit(1);
  }

  return snp_alleles;
}

Conf DetermineConf(std::vector<std::string> const &snp_seqs,
                   size_t site0,
                   size_t site1) {
  SnpAlleles snp_alleles0 = DetermineAlleles(snp_seqs, site0);
  SnpAlleles snp_alleles1 = DetermineAlleles(snp_seqs, site1);

  uint8_t a = snp_alleles0[0];
  uint8_t A = snp_alleles0[1];
  uint8_t b = snp_alleles1[0];
  uint8_t B = snp_alleles1[1];

  Conf conf;
  for (size_t hap = 0; hap < snp_seqs.size(); ++hap) {
    uint8_t base0 = ConvertCharToUInt8(snp_seqs[hap][site0]);
    uint8_t base1 = ConvertCharToUInt8(snp_seqs[hap][site1]);
    if (base0 == a && base1 == kUndefAllele) {
      ++conf.a_;
    } else if (base0 == A && base1 == kUndefAllele) {
      ++conf.A_;
    } else if (base0 == kUndefAllele && base1 == b) {
      ++conf.b_;
    } else if (base0 == kUndefAllele && base1 == B) {
      ++conf.B_;
    } else if (base0 == a && base1 == b) {
      ++conf.ab_;
    } else if (base0 == a && base1 == B) {
      ++conf.aB_;
    } else if (base0 == A && base1 == b) {
      ++conf.Ab_;
    } else if (base0 == A && base1 == B) {
      ++conf.AB_;
    } else if (base0 == kUndefAllele && base1 == kUndefAllele) {
      // 2-locus hap has N for both loci, which is valid.
    } else {
      fprintf(stderr, "SNPs are not diallelic.\n");
      std::exit(1);
    }
  }

  return conf;
}

void CleanSnpAndPos(std::vector<std::string> *snp_seqs,
                    std::vector<uint64_t> *snp_pos) {
  assert(snp_seqs->size() > 0);
  assert((*snp_seqs)[0].size() > 0);

  assert(snp_pos->size() > 0);
  assert((*snp_seqs)[0].size() == snp_pos->size());

  std::vector<std::string::iterator> snp_iters(snp_seqs->size());
  for (size_t i = 0; i < snp_seqs->size(); ++i) {
    snp_iters[i] = (*snp_seqs)[i].begin();
  }

  std::vector<uint64_t>::iterator pos_iter = snp_pos->begin();

  while (pos_iter != snp_pos->end()) {
    std::vector<std::string::const_iterator>
      const_snp_iters(snp_iters.size());
    for (size_t i = 0; i < snp_iters.size(); ++i) {
      const_snp_iters[i] = snp_iters[i];
    }

    SnpAlleles snp_alleles = DetermineAllelesCareful(const_snp_iters);
    bool remove_site = false;
    for (size_t j = 0; j < snp_alleles.size(); ++j) {
      if (!ValidSnpAlleleP(snp_alleles[j])) {
        remove_site = true;
      }
    }
    if (remove_site) {
      for (size_t i = 0; i < snp_iters.size(); ++i) {
        snp_iters[i] = (*snp_seqs)[i].erase(snp_iters[i]);
      }
      pos_iter = snp_pos->erase(pos_iter);
    } else {
      for (size_t i = 0; i < snp_iters.size(); ++i) {
        ++snp_iters[i];
      }
      ++pos_iter;
    }
  }
}
