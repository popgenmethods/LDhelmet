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

#ifndef LDHELMET_COMMON_RHO_FINDER_H_
#define LDHELMET_COMMON_RHO_FINDER_H_

#include <vector>

#include <boost/tuple/tuple.hpp>

typedef boost::tuple<std::vector<double>, std::vector<double> > RhoGrid;

std::vector<double> GetRhoList(RhoGrid const &rho_grid);

RhoGrid ParseRhoRange(std::vector<double> const &rho_range);

class RhoFinder {
 public:
  RhoFinder();

  explicit RhoFinder(RhoGrid const& rhoGrid_);

  // Returns index of greatest rho value in the table less than or
  // equal to given rho.
  size_t GetRhoID(double rho_value) const;

  std::vector<double> rho_list_;

 private:
  RhoGrid rho_grid_;

  // Records where each change of resolution occurs.
  std::vector<size_t> index_splits_;
};

#endif  // LDHELMET_COMMON_RHO_FINDER_H_
