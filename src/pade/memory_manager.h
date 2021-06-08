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

#ifndef LDHELMET_PADE_MEMORY_MANAGER_H_
#define LDHELMET_PADE_MEMORY_MANAGER_H_

#include <stdint.h>

#include <cstdlib>
#include <string>
#include <vector>

#include "pade/subtable.h"
#include "pade/subtable-inl.h"

template<typename ConfGenType>
class MemoryManager {
 public:
  MemoryManager(Table *table, uint32_t max_degree);
  ~MemoryManager();

  void FillTableFromFile(uint32_t m, uint32_t u,
                         ConfGenType const &in_conf_gen);

  void WriteElementsToFile(uint32_t m, uint32_t u,
                           ConfGenType const& in_conf_gen);

  void LoadAllElementsNeeded(uint32_t m,
                             uint32_t u,
                             ConfGenType const &in_conf_gen);

  inline void LoadElementsIncrementally(uint32_t m,
                                        uint32_t u,
                                        ConfGenType const &in_conf_gen) {
    if (u >= 3 && m + 1 <= max_degree_ / 2) {
      FillTableFromFile(m + 1, u - 3, in_conf_gen);
    }
  }

  inline void FreeMU(uint32_t m, uint32_t u) {
    Vec7().swap((*table_)[m][u]);
  }

  inline void FreeElementsIncrementally(uint32_t m, uint32_t u) {
    if (m >= 2) {
      FreeMU(m - 2, u);
    }
  }

  Table *table_;
  uint32_t const max_degree_;
  std::string tmp_dir_;
  std::vector<std::string> file_list_;
};

#endif  // LDHELMET_PADE_MEMORY_MANAGER_H_
