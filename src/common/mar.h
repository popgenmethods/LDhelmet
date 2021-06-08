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

#ifndef LDHELMET_COMMON_MAR_H_
#define LDHELMET_COMMON_MAR_H_

#include <stdint.h>

#include <string>

#include "common/conf.h"

// Marginal allele counts.
class Mar {
 public:
  Mar();

  Mar(uint32_t aMar, uint32_t AMar, uint32_t bMar, uint32_t BMar);

  std::string GetString() const;

  inline uint32_t ComputeDegree() const {
    return a_mar_ + A_mar_ + b_mar_ + B_mar_;
  }

  inline bool operator==(Mar const &other) const {
    return a_mar_ == other.a_mar_ &&
           A_mar_ == other.A_mar_ &&
           b_mar_ == other.b_mar_ &&
           B_mar_ == other.B_mar_;
  }

  inline bool operator!=(Mar const &other) const {
    return !(*this == other);
  }

  uint32_t a_mar_, A_mar_, b_mar_, B_mar_;
};

#endif  // LDHELMET_COMMON_MAR_H_
