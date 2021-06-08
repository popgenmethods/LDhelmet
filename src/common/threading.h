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

#ifndef LDHELMET_COMMON_THREADING_H_
#define LDHELMET_COMMON_THREADING_H_

#include <vector>

#include <boost/thread.hpp>

template<typename GenType>
class TaskBlocks {
 public:
  TaskBlocks() { }

  TaskBlocks(std::vector<GenType> const &blocks,
             size_t num_elems)
      : blocks_(blocks), num_elems_(num_elems) { }

  GenType& operator[](size_t index) {
    return blocks_[index];
  }

  GenType const& operator[](size_t index) const {
    return blocks_[index];
  }

  size_t size() const {
    return blocks_.size();
  }

  std::vector<GenType> blocks_;
  size_t num_elems_;
};

// Partition given generator into blocks.
template<typename GenType>
TaskBlocks<GenType> PartitionTask(size_t num_partitions,
                                  GenType const &in_gen) {
  TaskBlocks<GenType> task_blocks(
      std::vector<GenType>(num_partitions + 1), 0);

  size_t &num_elems = task_blocks.num_elems_;
  {
    GenType gen = in_gen;
    while (!gen.end()) {
      ++gen;
      ++num_elems;
    }
  }

  std::vector<GenType> &blocks = task_blocks.blocks_;
  {
    size_t block_size = num_elems/num_partitions;
    GenType gen = in_gen;
    size_t index = 0;
    for (size_t partition = 0; partition < num_partitions; ++partition) {
      while (index < partition * block_size) {
        ++gen;
        ++index;
      }
      blocks[partition] = gen;
    }
    while (!gen.end()) {
      ++gen;
    }
    blocks[num_partitions] = gen;
  }

  return task_blocks;
}

// Multithread tasks.
template<typename TaskMasterType, typename GenType>
size_t RunThreaded(size_t num_threads,
                   TaskMasterType const &task_master,
                   TaskBlocks<GenType> const &task_blocks) {
    if (num_threads > 256) {
      fprintf(stderr,
              "Number of threads requested is greater than 256. "
              "There's probably a bug somewhere or the user requested "
              "a high number of threads. This safety check will have to "
              "be removed to use more than 256 threads.\n");
      std::exit(1);
    }

    if (num_threads == 0) {
      fprintf(stderr,
              "Number of threads requested is 0. There is a bug somewhere.\n");
      std::exit(1);
    }

    assert(task_blocks.size() == num_threads + 1);

    boost::thread_group tgroup;
    for (size_t thread = 0; thread < task_blocks.size() - 1; ++thread) {
        tgroup.create_thread(
            task_master(task_blocks[thread], task_blocks[thread + 1]));
    }
    tgroup.join_all();

    return task_blocks.num_elems_;
}

#endif  // LDHELMET_COMMON_THREADING_H_
