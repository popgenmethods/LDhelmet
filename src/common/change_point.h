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

#ifndef LDHELMET_COMMON_CHANGE_POINT_H_
#define LDHELMET_COMMON_CHANGE_POINT_H_

#include <string>

#include <boost/lexical_cast.hpp>

class ChangePoint {
 public:
  ChangePoint(size_t snp_id, double rate) : snp_id_(snp_id), rate_(rate) { }

  inline
  std::string GetString() const {
    return '['
         + boost::lexical_cast<std::string>(snp_id_)
         + ':'
         + boost::lexical_cast<std::string>(rate_)
         + ']';
  }

  bool operator==(ChangePoint const &other) const {
    return snp_id_ == other.snp_id_ && rate_ == other.rate_;
  }

  size_t snp_id_;
  double rate_;
};

#endif  // LDHELMET_COMMON_CHANGE_POINT_H_
