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

#include "rjmcmc/acceptance_log.h"

#include <stdio.h>

AcceptanceLog AcceptanceLog::MergeStats(AcceptanceLog const &other) const {
  AcceptanceLog result;

  result.num_change_accepted_ = num_change_accepted_
                              + other.num_change_accepted_;
  result.num_change_proposed_ = num_change_proposed_
                              + other.num_change_proposed_;
  result.num_extend_accepted_ = num_extend_accepted_
                              + other.num_extend_accepted_;
  result.num_extend_proposed_ = num_extend_proposed_
                              + other.num_extend_proposed_;
  result.num_split_accepted_ = num_split_accepted_
                              + other.num_split_accepted_;
  result.num_split_proposed_ = num_split_proposed_
                              + other.num_split_proposed_;
  result.num_merge_accepted_ = num_merge_accepted_
                              + other.num_merge_accepted_;
  result.num_merge_proposed_ = num_merge_proposed_
                              + other.num_merge_proposed_;

  return result;
}

void AcceptanceLog::PrintStats() const {
  printf("accepted rate changes: %d/%d (%.2f)\n",
         num_change_accepted_, num_change_proposed_, FracChangeAccepted());
  printf("accepted extends: %d/%d (%.2f)\n",
         num_extend_accepted_, num_extend_proposed_, FracExtendAccepted());
  printf("accepted splits: %d/%d (%.2f)\n",
         num_split_proposed_, num_split_proposed_, FracSplitAccepted());
  printf("accepted merges: %d/%d (%.2f)\n",
         num_merge_accepted_, num_merge_proposed_, FracMergeAccepted());
  printf("--------\n");
  printf("total accepted: %d/%d (%.2f)\n",
         TotalAccepted(), TotalProposed(), FracAccepted());
}
