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

#ifndef LDHELMET_COMMON_CONF_GEN_H_
#define LDHELMET_COMMON_CONF_GEN_H_

#include <stdint.h>

#include <algorithm>
#include <cassert>

#include "common/conf.h"
#include "common/mar.h"
#include "common/predicates.h"

// Base class for "marinal allele counts" generators.
class MarGenBase {
 public:
  typedef Mar Elem;

  explicit MarGenBase(uint8_t degree);

  inline bool end() const {
    return end_;
  }

  inline Mar &operator*() {
    return mar_;
  }

  inline Mar const &operator*() const {
    return mar_;
  }

  inline bool operator==(MarGenBase const &other) const {
    return end_       == other.end_     &&
           degree_    == other.degree_  &&
           max_a_mar_ == other.max_a_mar_ &&
           max_A_mar_ == other.max_A_mar_ &&
           max_b_mar_ == other.max_b_mar_ &&
           max_B_mar_ == other.max_B_mar_ &&
           mar_       == other.mar_;
  }

  inline bool operator!=(MarGenBase const &other) const {
    return !(*this == other);
  }

 protected:
  bool end_;

  uint8_t degree_;
  Mar mar_;
  uint8_t max_a_mar_, max_A_mar_, max_b_mar_, max_B_mar_;
};

// Base case 1:
// Configurations where one ancestral allele remains at locus 1.
// aMar == 1 && AMar == 0.
class MarGenBaseA : public MarGenBase {
 public:
  MarGenBaseA();
  explicit MarGenBaseA(uint8_t degree);
  MarGenBaseA &operator++();
};

// Base case 2:
// Configurations where one ancestral allele remains at locus 2.
// bMar == 1 && BMar == 0.
class MarGenBaseB : public MarGenBase {
 public:
  MarGenBaseB();
  explicit MarGenBaseB(uint8_t degree);
  MarGenBaseB &operator++();
};

// Configurations that allow 1 mutation at locus 1 and 0 mutations at locus 2
// Iterates over all configurations where
//   0 <= aMar <= degree - 1
//   AMar == 1
//   0 <= bMar <= degree - 1 - aMar
//   BMar == degree - aMar - AMar - bMar (0 <= BMar <= degree - 1)
// Generate marginals for mutation at locus 1.
class MarGenMutA : public MarGenBase {
 public:
  MarGenMutA();
  explicit MarGenMutA(uint8_t degree);
  MarGenMutA &operator++();
};

// Configurations that allow 1 mutation at locus 2 and 0 mutations at locus 1
// Iterates over all configurations where
//   0 <= aMar <= degree - 1
//   0 <= AMar <= degree - 1 - aMar
//   bMar == degree - 1 - aMar - AMar (0 <= bMar <= degree - 1)
//   BMar == 1
// Generate marginals for mutation at locus 2.
class MarGenMutB : public MarGenBase {
 public:
  MarGenMutB();
  explicit MarGenMutB(uint8_t degree);
  MarGenMutB &operator++();
};

// Configurations that allow mutations at both loci
// Iterates over all configurations where
//   0 <= aMar <= degree - 2
//   AMar == 1
//   bMar == degree - 2 - aMar
//   BMar == 1
// Generate marginals for mutation at both loci.
class MarGenMutAB : public MarGenBase {
 public:
  MarGenMutAB();
  explicit MarGenMutAB(uint8_t degree);
  MarGenMutAB &operator++();
};

// Non-base case configurations that allow 0 mutations.
// Generate all marginals of a given degree.
class MarGen : public MarGenBase {
 public:
  MarGen();
  explicit MarGen(uint8_t degree);
  MarGen &operator++();
};

// Filter elements using a predicate, given a generator.
template<typename GenType,
         typename PredType = Otherwise<typename GenType::Elem> >
class GenPred {
 public:
  typedef typename GenType::Elem Elem;
  typedef GenType NestedGen;
  typedef PredType Pred;

  GenPred(GenType const &inGen = GenType(),
          PredType const &pred = PredType());

  GenPred &operator++();

  inline Elem &operator*() {
    return elem_;
  }

  inline Elem const &operator*() const {
    return elem_;
  }

  // Predicate member is not compared in equality comparison.
  inline bool operator==(GenPred const &other) const {
    return end_  == other.end_ &&
           gen_  == other.gen_ &&
           elem_ == other.elem_;
  }

  inline bool operator!=(GenPred const &other) const {
    return !(*this == other);
  }

  inline bool end() const {
    return end_;
  }

  inline GenType const &GetNestedGen() const {
    return gen_;
  }

  inline PredType const &GetPredicate() const {
    return pred_;
  }

 private:
  bool end_;

  GenType gen_;    // Nested generator.
  PredType pred_;  // Predicate.
  Elem elem_;      // Current element.
};

template<typename GenType, typename PredType>
GenPred<GenType, PredType>::GenPred(GenType const &gen,
                                    PredType const &pred)
    : end_(false), gen_(gen), pred_(pred) {
  while (!gen_.end()) {
    elem_ = *gen_;
    if (pred_(elem_)) {
      break;
    }
    ++gen_;
  }
  if (gen_.end()) {
    end_ = true;
  }
}

template<typename GenType, typename PredType>
GenPred<GenType, PredType> &GenPred<GenType, PredType>::operator++() {
  assert(!end());
  assert(!gen_.end());

  if (end()) {
    return *this;
  }

  do {
    ++gen_;
    elem_ = *gen_;
  } while (!gen_.end() && !pred_(elem_));

  if (gen_.end()) {
    end_ = true;
  }

  return *this;
}

// Generate all configurations, from degree 1 to degree max_degree given m.
// Used to compute g table for Pade coefficients.
class ConfGenM {
 public:
  typedef Conf Elem;

  ConfGenM();

  ConfGenM(uint8_t max_degree, uint8_t m);

  void Init(uint8_t max_degree, uint8_t m);

  ConfGenM &operator++();

  inline Conf &operator*() {
    return conf_;
  }

  inline Conf const &operator*() const {
    return conf_;
  }

  inline bool operator==(ConfGenM const &other) const {
    return conf_        == other.conf_       &&
           end_         == other.end_        &&
           max_degree_  == other.max_degree_ &&
           m_           == other.m_          &&
           degree_      == other.degree_     &&
           max_aB_      == other.max_aB_     &&
           max_Ab_      == other.max_Ab_     &&
           max_a_       == other.max_a_      &&
           max_A_       == other.max_A_      &&
           max_b_       == other.max_b_;
  }

  inline bool operator!=(ConfGenM const &other) const {
    return !(*this == other);
  }

  inline bool end() const {
    return end_;
  }

 private:
  bool end_;

  uint8_t max_degree_;
  uint8_t m_;

  uint8_t degree_;
  Conf conf_;

  uint8_t max_aB_, max_Ab_, max_a_, max_A_, max_b_;
};

// Generate all subconfigurations of c given m.
// Used to compute q values for Pade coefficients.
class ConfGenMC {
 public:
  ConfGenMC(Conf const &c, uint8_t m);
  ConfGenMC &operator++();

  inline Conf &operator*() {
      assert(conf_.ComputeDegree() == degree_);
      assert(conf_.ComputeDegree() >= (uint32_t)(2 * m_)
          && conf_.ComputeDegree() >= 1);
      return conf_;
  }

  inline Conf const &operator*() const {
      return conf_;
  }

  inline bool end() const {
      return end_;
  }

  inline bool operator==(ConfGenMC const &other) const {
      return conf_   == other.conf_   &&
             end_    == other.end_    &&
             c_      == other.c_      &&
             m_      == other.m_      &&
             degree_ == other.degree_ &&
             max_ab_  == other.max_ab_  &&
             max_aB_  == other.max_aB_  &&
             max_Ab_  == other.max_Ab_;
  }

  inline bool operator!=(ConfGenMC const &other) const {
      return !(*this == other);
  }

 private:
  bool end_;

  Conf const c_;
  uint8_t m_;

  uint8_t degree_;
  Conf conf_;

  uint8_t max_ab_, max_aB_, max_Ab_;
};

// This generator is inclusive of the start ID and exclusive of the end ID.
class ConfIDGen {
 public:
  typedef uint64_t Elem;

  ConfIDGen();

  ConfIDGen(uint64_t start_id, uint64_t end_id);

  inline bool end() const {
    return end_;
  }

  ConfIDGen &operator++();

  inline uint64_t &operator*() {
    return cur_conf_id_;
  }

  inline uint64_t const &operator*() const {
    return cur_conf_id_;
  }

  inline bool operator==(ConfIDGen const &other) const {
    return cur_conf_id_ == other.cur_conf_id_ &&
           end_         == other.end_         &&
           start_id_    == other.start_id_    &&
           end_id_      == other.end_id_;
  }

  inline bool operator!=(ConfIDGen const &other) const {
    return !(*this == other);
  }

 private:
  bool end_;
  uint64_t start_id_, end_id_;
  uint64_t cur_conf_id_;
};

#endif  // LDHELMET_COMMON_CONF_GEN_H_
