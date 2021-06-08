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

#ifndef LDHELMET_PADE_SUBTABLE_H_
#define LDHELMET_PADE_SUBTABLE_H_

#include <stddef.h>
#include <stdint.h>

#include <vector>

#include "common/conf_gen.h"
#include "common/vector_definitions.h"

template<typename T> class Subtable;

typedef Subtable<Vec7> Table;

template<typename T>
class Subtable {
 public:
  Subtable() { }

  Subtable(size_t size1, size_t size2);

  Subtable(size_t size1, size_t size2, T init_val);

  inline std::vector<T>& operator[](uint32_t m) {
    return subtable_[m/2];
  }

  inline std::vector<T> const& operator[](uint32_t m) const {
    return subtable_[m/2];
  }

  inline std::vector<std::vector<T> > &GetSubtable() {
    return subtable_;
  }

  inline std::vector<std::vector<T> > const &GetSubtable() const {
    return subtable_;
  }

  // ab, aB, Ab, a, A, b, B.
  // AB is inferred from m and ab, aB, Ab.
  void AllocateMemory(uint32_t max_degree, uint32_t m, uint32_t u);

 private:
  std::vector<std::vector<T> > subtable_;
};

double GetTable(Table const &table, uint32_t m, uint32_t u, Conf const &conf);

void SetTable(Table *table, uint32_t m, uint32_t u,
              Conf const &conf, double value);

#endif  // LDHELMET_PADE_SUBTABLE_H_
