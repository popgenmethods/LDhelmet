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

#ifndef LDHELMET_COMMON_BINARY_SEARCH_H_
#define LDHELMET_COMMON_BINARY_SEARCH_H_

#include <stddef.h>

#include <cassert>
#include <utility>
#include <vector>

// Comparison function object.
template<typename T>
struct DefaultCompare {
  bool operator()(T const &a, T const &b) const {
      return a < b;
  }
};

// Binary search.
// Assumes xs is in strictly increasing order.
template<typename T, typename Compare>
inline std::pair<size_t, size_t> BinarySearch(std::vector<T> const &xs,
                                              T const &x,
                                              Compare const &compare) {
  size_t left_pt = 0;
  size_t right_pt = xs.size();

  while (right_pt - left_pt > 1) {
    size_t mid_pt = (right_pt + left_pt)/2;
    if (compare(x, xs[mid_pt])) {
      right_pt = mid_pt;
    } else {
      left_pt = mid_pt;
    }
  }

  assert(right_pt - left_pt == 1);

  return std::pair<size_t, size_t>(left_pt, right_pt);
}

template<typename T>
inline std::pair<size_t, size_t> BinarySearch(std::vector<T> const &xs,
                                              T const &x) {
  return BinarySearch(xs, x, DefaultCompare<T>());
}

#endif  // LDHELMET_COMMON_BINARY_SEARCH_H_
