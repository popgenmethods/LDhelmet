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

#include "common/mar.h"

#include <string>

#include <boost/lexical_cast.hpp>

#include "common/conf.h"

Mar::Mar() : a_mar_(0), A_mar_(0), b_mar_(0), B_mar_(0) { }

Mar::Mar(uint32_t a_mar, uint32_t A_mar, uint32_t b_mar, uint32_t B_mar)
    : a_mar_(a_mar), A_mar_(A_mar), b_mar_(b_mar), B_mar_(B_mar) { }

std::string Mar::GetString() const {
  return boost::lexical_cast<std::string>(a_mar_) + ' ' +
         boost::lexical_cast<std::string>(A_mar_) + ' ' +
         boost::lexical_cast<std::string>(b_mar_) + ' ' +
         boost::lexical_cast<std::string>(B_mar_);
}
