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

#ifndef LDHELMET_PADE_MEMORY_MANAGER_INL_H_
#define LDHELMET_PADE_MEMORY_MANAGER_INL_H_

#include <stdint.h>
#include <stdio.h>

#include <string>

#include <boost/lexical_cast.hpp>

template<typename ConfGenType>
MemoryManager<ConfGenType>::MemoryManager(Table *table, uint32_t max_degree)
  : table_(table), max_degree_(max_degree) {
    // Create temporary directory.
    char tmp_dir_name_template[] = "tmp_pade_XXXXXX";
    tmp_dir_ = mkdtemp(tmp_dir_name_template);
    if (tmp_dir_ == "") {
      fprintf(stderr, "Unable to create temporary directory.\n");
      std::exit(1);
    }
  }

template<typename ConfGenType>
MemoryManager<ConfGenType>::~MemoryManager() {
  printf("Removing temporary files.\n");

  // Remove each temporary file.
  for (int i = 0; i < static_cast<int>(file_list_.size()); ++i) {
    int return_status = unlink(file_list_[i].c_str());
    if (return_status != 0) {
      fprintf(stderr,
              "Error deleting temporary file in MemoryManager's destructor.\n");
    }
  }

  // Remove the temporary directory.
  {
    int return_status = rmdir(tmp_dir_.c_str());
    if (return_status != 0) {
      fprintf(stderr,
              "Error deleting temporary directory in MemoryManager's "
              "destructor.\n");
    }
  }
}

template<typename ConfGenType>
void MemoryManager<ConfGenType>::FillTableFromFile(uint32_t m, uint32_t u,
    ConfGenType const &in_conf_gen) {
  // Construct tmp file name to read from.
  std::string tmp_file_in = "tmp_g_"
                          + boost::lexical_cast<std::string>(m)
                          + "_"
                          + boost::lexical_cast<std::string>(u);
  std::string tmp_path_in = tmp_dir_ + "/" + tmp_file_in;
  FILE *fp = fopen(tmp_path_in.c_str(), "rb");
  if (fp == NULL) {
    fprintf(stderr,
            "Could not open temporary file for reading: %s.\n",
            tmp_path_in.c_str());
    std::exit(1);
  }

  table_->AllocateMemory(max_degree_, m, u);

  // Change the conf_gen's max_degree_ and m_, but keep the same
  // predicate as the one passed in.
  ConfGenType conf_gen(typename ConfGenType::NestedGen(max_degree_, m),
      in_conf_gen.GetPredicate());
  while (!conf_gen.end()) {
    double g_value;

    int num_read = fread(reinterpret_cast<char *>(&g_value),
                         sizeof(g_value), 1, fp);
    assert(num_read == 1);
    SetTable(table_, m, u, *conf_gen, g_value);

    ++conf_gen;
  }
  assert(fp != NULL);
  fclose(fp);
}

template<typename ConfGenType>
void MemoryManager<ConfGenType>::WriteElementsToFile(
    uint32_t m, uint32_t u,
    ConfGenType const& in_conf_gen) {
  // Construct tmp file name to write to.
  std::string tmp_file_out = "tmp_g_"
                           + boost::lexical_cast<std::string>(m)
                           + "_"
                           + boost::lexical_cast<std::string>(u);
  std::string tmp_path_out = tmp_dir_ + "/" + tmp_file_out;
  FILE *fp = fopen(tmp_path_out.c_str(), "wb");
  if (fp == NULL) {
    fprintf(stderr,
            "Could not open temporary file for writing: %s.\n",
            tmp_path_out.c_str());
    std::exit(1);
  }

  file_list_.push_back(tmp_path_out);

  ConfGenType conf_gen = in_conf_gen;
  while (!conf_gen.end()) {
    Conf conf = *conf_gen;
    double g_value = GetTable(*table_, m, u, conf);
    int num_written = fwrite(reinterpret_cast<char *>(&g_value),
                             sizeof(g_value), 1, fp);
    assert(num_written == 1);
    ++conf_gen;
  }
  assert(fp != NULL);
  fclose(fp);
}

template<typename ConfGenType>
void MemoryManager<ConfGenType>::LoadAllElementsNeeded(
    uint32_t m,
    uint32_t u,
    ConfGenType const &in_conf_gen) {
  if (u >= 3 && m + 1 <= max_degree_ / 2) {
    FillTableFromFile(m + 1, u - 3, in_conf_gen);
  }
  if (u >= 1) {
    FillTableFromFile(m, u - 2, in_conf_gen);
  }
  if (m >= 1 && u >= 1) {
    FillTableFromFile(m - 1, u - 1, in_conf_gen);
  }
  if (m >= 2) {
    FillTableFromFile(m - 2, u, in_conf_gen);
  }
}

#endif  // LDHELMET_PADE_MEMORY_MANAGER_INL_H_
