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

#ifndef LDHELMET_PADE_SUBTABLE_INL_H_
#define LDHELMET_PADE_SUBTABLE_INL_H_

#include <algorithm>
#include <vector>

template<typename T>
Subtable<T>::Subtable(size_t size1,
                      size_t size2)
    : subtable_(size1/2 + 1, std::vector<T>(size2, T())) { }

template<typename T>
Subtable<T>::Subtable(size_t size1,
                      size_t size2,
                      T init_val)
    : subtable_(size1/2 + 1, std::vector<T>(size2, init_val)) { }

template<typename T>
void Subtable<T>::AllocateMemory(uint32_t max_degree, uint32_t m, uint32_t u) {
  uint64_t num_elem = 0;

  assert(max_degree >= 1);
  Subtable &s_table = *this;
  uint32_t degree_limit = max_degree
                        - std::max(static_cast<uint32_t>(1), 2*m) + 1;
  s_table[m][u] = Vec7(degree_limit);
  // degree is (true degree) - std::max(1, 2*m).
  // The true degree should range from 2*m to max_degree if m >= 1
  //                              from 1 to max_degree if m == 0.
  // The fake degree should range from 0 to max_degree - 2*m if m >= 1
  //                              from 0 to max_degree - 1 if m == 0.
  for (uint32_t degree = 0; degree < degree_limit; ++degree) {
    s_table[m][u][degree].resize(m+1);
    for (uint32_t ab = 0; ab <= m; ++ab) {
      uint32_t max_aB = m - ab;
      s_table[m][u][degree][ab].resize(max_aB + 1);
      for (uint32_t aB = 0; aB <= max_aB; ++aB) {
        uint32_t max_Ab = max_aB - aB;
        s_table[m][u][degree][ab][aB].resize(max_Ab + 1);
        for (uint32_t Ab = 0; Ab <= max_Ab; ++Ab) {
          // (fake degree) + std::max(1,2*m) == (true degree),
          //   so max_a = (true degree) - 2*m.
          uint32_t max_a = degree
                         + std::max(static_cast<uint32_t>(1), 2*m) - 2*m;
          s_table[m][u][degree][ab][aB][Ab].resize(max_a + 1);
          for (uint32_t a = 0; a <= max_a; ++a) {
            uint32_t max_A = max_a - a;
            s_table[m][u][degree][ab][aB][Ab][a].resize(max_A + 1);
            for (uint32_t A = 0; A <= max_A; ++A) {
              uint32_t max_b = max_A - A;
              s_table[m][u][degree][ab][aB][Ab][a][A].resize(max_b + 1);
              num_elem += max_b + 1;
            }
          }
        }
      }
    }
  }
}

#endif  // LDHELMET_PADE_SUBTABLE_INL_H_
