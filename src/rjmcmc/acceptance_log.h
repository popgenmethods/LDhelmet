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

#ifndef LDHELMET_RJMCMC_ACCEPTANCE_LOG_H_
#define LDHELMET_RJMCMC_ACCEPTANCE_LOG_H_

#include <stdint.h>

class AcceptanceLog {
 public:
  AcceptanceLog()
      : num_change_proposed_(),
        num_change_accepted_(),
        num_extend_proposed_(),
        num_extend_accepted_(),
        num_split_proposed_(),
        num_split_accepted_(),
        num_merge_proposed_(),
        num_merge_accepted_() { }

  inline void AcceptChange(bool accept_p) {
    if (accept_p) {
      ++num_change_accepted_;
    }
    ++num_change_proposed_;
  }

  inline void AcceptExtend(bool accept_p) {
    if (accept_p) {
      ++num_extend_accepted_;
    }
    ++num_extend_proposed_;
  }

  inline void AcceptSplit(bool accept_p) {
    if (accept_p) {
      ++num_split_accepted_;
    }
    ++num_split_proposed_;
  }

  inline void AcceptMerge(bool accept_p) {
    if (accept_p) {
      ++num_merge_accepted_;
    }
    ++num_merge_proposed_;
  }

  inline double FracChangeAccepted() const {
    return static_cast<double>(num_change_accepted_)
           / static_cast<double>(num_change_proposed_);
  }

  inline double FracExtendAccepted() const {
    return static_cast<double>(num_extend_accepted_)
           / static_cast<double>(num_extend_proposed_);
  }

  inline double FracSplitAccepted() const {
    return static_cast<double>(num_split_accepted_)
           / static_cast<double>(num_split_proposed_);
  }

  inline double FracMergeAccepted() const {
      return static_cast<double>(num_merge_accepted_)
             / static_cast<double>(num_merge_proposed_);
  }

  inline uint32_t TotalAccepted() const {
      return num_change_accepted_
           + num_extend_accepted_
           + num_split_accepted_
           + num_merge_accepted_;
  }

  inline uint32_t TotalProposed() const {
      return num_change_proposed_
           + num_extend_proposed_
           + num_split_proposed_
           + num_merge_proposed_;
  }

  inline double FracAccepted() const {
      return static_cast<double>(TotalAccepted())
             / static_cast<double>(TotalProposed());
  }

  // Combine stats from two AcceptanceLogs.
  AcceptanceLog MergeStats(AcceptanceLog const &other) const;

  void PrintStats() const;

  uint32_t num_change_proposed_;
  uint32_t num_change_accepted_;

  uint32_t num_extend_proposed_;
  uint32_t num_extend_accepted_;

  uint32_t num_split_proposed_;
  uint32_t num_split_accepted_;

  uint32_t num_merge_proposed_;
  uint32_t num_merge_accepted_;
};

#endif  // LDHELMET_RJMCMC_ACCEPTANCE_LOG_H_
