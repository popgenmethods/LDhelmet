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

#ifndef LDHELMET_FIND_CONFS_FIND_CONFS_H_
#define LDHELMET_FIND_CONFS_FIND_CONFS_H_

#include <set>
#include <string>
#include <vector>

#include <boost/thread.hpp>
#include <boost/tuple/tuple.hpp>

#include "common/conf.h"

boost::tuple<std::vector<std::string>, std::vector<uint64_t> > FindSnpsHelper(
    std::vector<std::string> const &seqs,
    size_t begin,
    size_t end);

boost::tuple<std::vector<std::string>, std::vector<uint64_t> > FindSnps(
    uint32_t const num_threads,
    std::vector<std::string> const &seqs);

std::set<Conf> FindConfs(uint32_t const window_size,
                         std::vector<std::string> const &snp_seqs);

std::set<Conf> MergeConfs(std::vector<std::set<Conf> > const &conf_sets);
class FindSnpsPartition {
 public:
  FindSnpsPartition(std::vector<std::string> const &seqs,
                    size_t begin,
                    size_t end,
                    boost::tuple<std::vector<std::string>,
                                 std::vector<uint64_t> > &result)
      : seqs_(seqs), begin_(begin), end_(end), result_(result) { }

  void operator()() {
    result_ = FindSnpsHelper(seqs_, begin_, end_);
  }

  std::vector<std::string> const &seqs_;
  size_t begin_, end_;
  boost::tuple<std::vector<std::string>, std::vector<uint64_t> > &result_;
};

struct FindConfsPartition {
  FindConfsPartition(uint32_t window_size,
                     std::vector<std::string> const &seqs,
                     std::set<Conf> &result)
      : window_size_(window_size), seqs_(seqs), result_(result) { }

  void operator()() {
    result_ = FindConfs(window_size_, seqs_);
  }

  uint32_t window_size_;
  std::vector<std::string> const &seqs_;
  std::set<Conf> &result_;
};

std::set<Conf> FindConfsMT(uint32_t const num_threads,
                           uint32_t const window_size,
                           std::vector<std::string> const &snp_seqs);

std::vector<Conf> FindConfsFromFiles(
    uint32_t const num_threads,
    uint32_t const window_size,
    std::vector<std::string> const &input_files);

#endif  // LDHELMET_FIND_CONFS_FIND_CONFS_H_
