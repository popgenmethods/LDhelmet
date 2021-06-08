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

#ifndef LDHELMET_COMMON_PREDICATES_H_
#define LDHELMET_COMMON_PREDICATES_H_

#include <stdint.h>

#include "common/conf.h"
#include "common/mar.h"

// Predicate to determine whether a configuration can exist in real data for a
// given sample size
// For example, with n=10, the configuration 20 0 0 0 0 0 0 0 cannot exist in
// real data because there can be at most 10 alleles defined for the first
// locus (since n=10) and at most 10 alleles for the second locus
// (configuration: a A b B ab aB Ab AB)
// Neither is this configuration valid: 10 10 0 0 0 0 0 0, for the same reason
// as above
// But this is valid: 10 0 10 0 0 0 0 0
// A configuration with only ancestral alleles is valid, e.g. 10 0 10 0 0 0 0 0
// But a configuration that is missing an ancestral allele for one of the loci
// is not valid, e.g. 0 10 10 0 0 0 0 0 is not valid because it has no ancestral
// allele for the first locus
class SCCMutation {
 public:
  SCCMutation() : max_locus_(-1) { }

  explicit SCCMutation(uint32_t max_locus) : max_locus_(max_locus) { }

  inline bool operator()(Mar const &mar) const {
    return mar.a_mar_ + mar.A_mar_ <= max_locus_ &&
           mar.b_mar_ + mar.B_mar_ <= max_locus_;
  }

  inline bool operator()(Conf const &conf) const {
    return operator()(Mar(conf.a_ + conf.ab_ + conf.aB_,
                          conf.A_ + conf.Ab_ + conf.AB_,
                          conf.b_ + conf.ab_ + conf.Ab_,
                          conf.B_ + conf.aB_ + conf.AB_));
  }

  uint32_t max_locus_;
};

class SCCAll {
 public:
  inline bool operator()(Mar const &mar) const {
    return true;
  }

  inline bool operator()(Conf const &conf) const {
    return true;
  }
};

template<typename ArgType>
class Otherwise {
 public:
  typedef ArgType Arg;

  inline bool operator()(ArgType const &x) const {
    return true;
  }
};

template<typename ArgType, typename PredType = Otherwise<ArgType> >
class PredBase {
 public:
  typedef ArgType Arg;
  typedef PredType Tail;

  PredBase() : pred_(PredType()) { }
  PredBase(PredType const &pred) : pred_(pred) { }

  PredType pred_;
};

template<typename PredType = Otherwise<Mar> >
class DecideToSolve : PredBase<Mar, PredType> {
 public:
  DecideToSolve()
      : PredBase<Mar, PredType>(PredType()) { }

  DecideToSolve(PredType const &pred)
      : PredBase<Mar, PredType>(pred) { }

  bool operator()(Mar const &mar) const {
    return (mar.a_mar_ > 0 && mar.b_mar_ > 0) &&
           (mar.a_mar_ == mar.b_mar_ ? mar.A_mar_ >= mar.B_mar_
                                     : mar.a_mar_ >  mar.b_mar_) &&
           // Evaluate other predicates.
           PredBase<Mar, PredType>::pred_(mar);
  }
};

template<typename PredType = Otherwise<Mar> >
class BaseCaseAPred : PredBase<Mar, PredType> {
 public:
  BaseCaseAPred()
      : PredBase<Mar, PredType>(PredType()) { }

  BaseCaseAPred(PredType const &pred)
      : PredBase<Mar, PredType>(pred) { }

  inline bool operator()(Mar const &mar) const {
    return PredBase<Mar, PredType>::pred_(mar);
  }
};

template<typename PredType = Otherwise<Mar> >
class BaseCaseBPred : PredBase<Mar, PredType> {
 public:
  BaseCaseBPred()
      : PredBase<Mar, PredType>(PredType()) { }

  BaseCaseBPred(PredType const &pred)
      : PredBase<Mar, PredType>(pred) { }

  inline bool operator()(Mar const &mar) const {
           // Skip locus 1 base case configurations.
    return !(mar.a_mar_ == 1 && mar.A_mar_ == 0) &&
           // Evaluate other predicates.
           PredBase<Mar, PredType>::pred_(mar);
  }
};

template<typename PredType = Otherwise<Mar> >
class NoMutPred : PredBase<Mar, PredType> {
 public:
  NoMutPred()
      : PredBase<Mar, PredType>(PredType()) { }

  NoMutPred(PredType const &pred)
      : PredBase<Mar, PredType>(pred) { }

  inline bool operator()(Mar const &mar) const {
           // Skip locus 1 base case configurations.
    return !(mar.a_mar_ == 1 && mar.A_mar_ == 0) &&
           // Skip locus 2 base case configurations.
           !(mar.b_mar_ == 1 && mar.B_mar_ == 0) &&
           // Skip mutation configurations.
           !(mar.A_mar_ == 1 || mar.B_mar_ == 1) &&
           // Evaluate other predicates.
           PredBase<Mar, PredType>::pred_(mar);
  }
};

template<typename PredType = Otherwise<Mar> >
class MutAPred : PredBase<Mar, PredType> {
 public:
  MutAPred()
      : PredBase<Mar, PredType>(PredType()) { }

  MutAPred(PredType const &pred)
      : PredBase<Mar, PredType>(pred) { }

  inline bool operator()(Mar const &mar) const {
           // Skip locus 1 base case configurations.
    return !(mar.a_mar_ == 1 && mar.A_mar_ == 0) &&
           // Skip locus 2 base case configurations.
           !(mar.b_mar_ == 1 && mar.B_mar_ == 0) &&
           // Skip two mutation configurations.
           !(mar.A_mar_ == 1 && mar.B_mar_ == 1) &&
           // Evaluate other predicates.
           PredBase<Mar, PredType>::pred_(mar);
  }
};

template<typename PredType = Otherwise<Mar> >
class MutBPred : PredBase<Mar, PredType> {
 public:
  MutBPred()
      : PredBase<Mar, PredType>(PredType()) { }

  MutBPred(PredType const &pred)
      : PredBase<Mar, PredType>(pred) { }

  inline bool operator()(Mar const &mar) const {
           // Skip locus 1 base case configurations.
    return !(mar.a_mar_ == 1 && mar.A_mar_ == 0) &&
           // Skip locus 2 base case configurations.
           !(mar.b_mar_ == 1 && mar.B_mar_ == 0) &&
           // Skip two mutation configurations.
           !(mar.A_mar_ == 1 && mar.B_mar_ == 1) &&
           // Evaluate other predicates.
           PredBase<Mar, PredType>::pred_(mar);
  }
};

template<typename PredType = Otherwise<Mar> >
class MutABPred : PredBase<Mar, PredType> {
 public:
  MutABPred()
      : PredBase<Mar, PredType>(PredType()) { }

  MutABPred(PredType const &pred)
      : PredBase<Mar, PredType>(pred) { }

  inline bool operator()(Mar const &mar) const {
    return PredBase<Mar, PredType>::pred_(mar);
  }
};

#endif  // LDHELMET_COMMON_PREDICATES_H_
