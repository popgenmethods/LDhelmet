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

#ifndef LDHELMET_COMMON_CONF_H
#define LDHELMET_COMMON_CONF_H

#include <stddef.h>
#include <stdint.h>

#include <functional>
#include <string>

#include <boost/functional/hash.hpp>

// haplotype configuration
class Conf {
 public:
  class BinaryRep {
   public:
    typedef uint64_t RepType;

    BinaryRep();
    explicit BinaryRep(Conf const &conf);

    inline RepType rep() const {
      return rep_;
    }

    RepType rep_;
  };

  Conf();

  Conf(uint8_t a, uint8_t A, uint8_t b, uint8_t B,
       uint8_t ab, uint8_t aB, uint8_t Ab, uint8_t AB);

  explicit Conf(BinaryRep const &binaryRep);

  std::string GetString() const;

  Conf Adda(int i, int quantity) const;
  Conf Addb(int j, int quantity) const;
  Conf Addr(int i, int j, int quantity) const;

  bool operator<(Conf const &conf2) const;

  inline int Geta(int i) const {
    if (i == 0) {
      return a_;
    } else {
      return A_;
    }
  }

  inline int Getb(int j) const {
    if (j == 0) {
      return b_;
    } else {
      return B_;
    }
  }

  inline int Getr(int i, int j) const {
    if (i == 0) {
      if (j == 0) {
        return ab_;
      } else {
        return aB_;
      }
    } else {
      if (j == 0) {
        return Ab_;
      } else {
        return AB_;
      }
    }
  }

  inline uint32_t ComputeaMar() const {
    return static_cast<uint32_t>(a_) + ab_ + aB_;
  }

  inline uint32_t ComputeAMar() const {
    return static_cast<uint32_t>(A_) + Ab_ + AB_;
  }

  inline uint32_t ComputebMar() const {
    return static_cast<uint32_t>(b_) + ab_ + Ab_;
  }

  inline uint32_t ComputeBMar() const {
    return static_cast<uint32_t>(B_) + aB_ + AB_;
  }

  inline uint32_t ComputeM() const {
    return static_cast<uint32_t>(ab_) + aB_ + Ab_ + AB_;
  }

  inline uint32_t ComputeDegree() const {
    return static_cast<uint32_t>(a_) + A_ + b_ + B_
         + 2*(static_cast<uint32_t>(ab_) + aB_ + Ab_ + AB_);
  }

  inline uint32_t ComputeSize() const {
    return static_cast<uint32_t>(a_) + A_ + b_ + B_
         + (static_cast<uint32_t>(ab_) + aB_ + Ab_ + AB_);
  }

  inline bool operator==(Conf const &other) const {
    return a_  == other.a_  &&
           A_  == other.A_  &&
           b_  == other.b_  &&
           B_  == other.B_  &&
           ab_ == other.ab_ &&
           aB_ == other.aB_ &&
           Ab_ == other.Ab_ &&
           AB_ == other.AB_;
  }

  inline bool operator!=(Conf const &other) const {
    return !(*this == other);
  }

  uint8_t a_, A_, b_, B_, ab_, aB_, Ab_, AB_;
};

struct ConfHash : std::unary_function<Conf, size_t> {
  inline size_t operator()(Conf const &conf) const {
    size_t seed = 0;
    boost::hash_combine(seed, conf.a_);
    boost::hash_combine(seed, conf.A_);
    boost::hash_combine(seed, conf.b_);
    boost::hash_combine(seed, conf.B_);
    boost::hash_combine(seed, conf.ab_);
    boost::hash_combine(seed, conf.aB_);
    boost::hash_combine(seed, conf.Ab_);
    boost::hash_combine(seed, conf.AB_);
    return seed;
  }
};

#endif  // LDHELMET_COMMON_CONF_H
