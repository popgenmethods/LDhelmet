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

#include <stdio.h>

#include <cassert>
#include <cmath>

#include "common/rho_finder.h"

RhoGrid ParseRhoRange(std::vector<double> const &rho_range) {
  if (rho_range.size() % 2 == 0 || rho_range.size() == 0) {
    fprintf(stderr,
            "Rho range is invalid.\n");
    exit(1);
  }

  RhoGrid rho_grid;
  std::vector<double> &rho_points = rho_grid.get<0>();
  std::vector<double> &deltas = rho_grid.get<1>();
  for (size_t i = 0; i + 1 < rho_range.size(); i += 2) {
    rho_points.push_back(rho_range[i]);
    deltas.push_back(rho_range[i + 1]);
  }
  rho_points.push_back(rho_range.back());

  return rho_grid;
}

RhoFinder::RhoFinder() { }

RhoFinder::RhoFinder(RhoGrid const& rho_grid)
    : rho_list_(GetRhoList(rho_grid)), rho_grid_(rho_grid) {
  std::vector<double> &rho_points = rho_grid_.get<0>();
  size_t rho_points_index = 0;
  for (size_t i = 0; i < rho_list_.size(); ++i) {
    if (rho_list_[i] == rho_points[rho_points_index]) {
      index_splits_.push_back(i);
      ++rho_points_index;
    }
  }
  assert(rho_points_index == rho_points.size());
  assert(index_splits_.size() == rho_points.size());
}

//returns the index of the largest rho in rho_grid
//that is less than or equal to rho_value
size_t RhoFinder::GetRhoID(double rho_value) const {
  assert(rho_value >= 0);

  std::vector<double> const &rho_points = rho_grid_.get<0>();
  std::vector<double> const &deltas = rho_grid_.get<1>();

  if (rho_value < rho_points[0]) {
    fprintf(stderr,
            "The request rho is less than all rho values in "
            "the lookup table. "
            "This indicates a problem with the lookup table. "
            "The lookup table must include rho value 0.0.\n");
    exit(1);
  }
    
  if(rho_value >= rho_list_[rho_list_.size() - 1] - (1e-10)*deltas[rho_points.size() - 2]){
     return rho_list_.size() - 1;
  }


  for (size_t i = 0; i < rho_points.size() - 1; ++i) {
    double truncated_rho_value = std::floor(rho_value / deltas[i]) * deltas[i];
      if (truncated_rho_value > rho_value) {
          //assert(truncated_rho_value - rho_value < (1e-10)*deltas[i]);
          truncated_rho_value = rho_value;
      }
      if (rho_value - truncated_rho_value > (1-1e-10)*deltas[i]){
          truncated_rho_value += deltas[i];
      }
      if (truncated_rho_value < rho_points[i+1]){
          int offset = static_cast<int>((truncated_rho_value - rho_points[i])
                                        / deltas[i] + 0.5);
          size_t rho_id = index_splits_[i] + offset;
          //assert(rho_id >= 0);
          //assert(rho_id < rho_list_.size());
          //assert(rho_value > rho_list_[rho_id] - (1e-10)*deltas[i]);
          //assert(rho_value < rho_list_[rho_id+1]);
          return rho_id;
      }
  }
  assert(false);        //should not get here.
  return rho_list_.size() - 1;
}

std::vector<double> GetRhoList(RhoGrid const &rho_grid) {
  assert(rho_grid.get<0>().size() > 0 && rho_grid.get<1>().size() >0);
  assert(rho_grid.get<0>().size() == rho_grid.get<1>().size() + 1);

  std::vector<double> const &rho_points = rho_grid.get<0>();
  std::vector<double> const &deltas = rho_grid.get<1>();

  std::vector<double> rho_list;
  for (size_t i = 0; i < rho_points.size() - 1; ++i) {
    double base_rho = rho_points[i];
    double delta = deltas[i];
    if (delta <= 0.0) {
      fprintf(stderr,
              "Rho range is invalid. Delta between rhos is non-positive.\n");
      exit(1);
    }
    if (rho_points[i + 1] <= rho_points[i]) {
      fprintf(stderr,
              "Rho range is invalid. "
              "Rho values are not monotonically increasing.\n");
      exit(1);
    }
    size_t rho_multiplier = 0;
    double cur_rho = base_rho + rho_multiplier * delta;
    while (cur_rho < rho_points[i + 1] - (1e-10)*delta) {
      rho_list.push_back(cur_rho);
      ++rho_multiplier;
      cur_rho = base_rho + static_cast<double>(rho_multiplier) * delta;
    }
  }
  rho_list.push_back(rho_points.back());

  assert(rho_list.size()>0);
  assert(rho_list[0] == rho_points[0]);

  for (size_t rho_id = 1; rho_id < rho_list.size(); ++rho_id) {
    if (rho_list[rho_id] < rho_list[rho_id - 1]) {
      fprintf(stderr,
              "Rho values must be in increasing order.");
      exit(1);
    }
  }

  return rho_list;
}
