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

#ifndef LDHELMET_PADE_G_TASK_H_
#define LDHELMET_PADE_G_TASK_H_

#include <stdint.h>

#include <vector>

#include "pade/compute_g.h"

template<typename ConfGenType>
class GTaskWorker {
 public:
  GTaskWorker(double theta,
              Table *table,
              uint32_t m,
              uint32_t u,
              ConfGenType const &start_conf_gen,
              ConfGenType const &end_conf_gen)
      : theta_(theta),
        table_(table),
        m_(m),
        u_(u),
        start_conf_gen_(start_conf_gen),
        end_conf_gen_(end_conf_gen) { }

  void operator()() {
    for (ConfGenType conf_gen = start_conf_gen_;
         conf_gen != end_conf_gen_;
         ++conf_gen) {
      assert(!conf_gen.end());
      ComputeG(theta_, table_, m_, u_, *conf_gen);
    }
  }

 private:
  double theta_;
  Table *table_;
  uint32_t m_, u_;
  ConfGenType start_conf_gen_, end_conf_gen_;
};

// Contains parameters for work; produces workers when given start and end conf
// iterators.
class GTaskMaster {
 public:
  GTaskMaster(double theta, Table *table, uint32_t m, uint32_t u)
      : theta_(theta),
        table_(table),
        m_(m),
        u_(u) { }

  template<typename ConfGenType>
  GTaskWorker<ConfGenType> operator()(ConfGenType const &start_conf_gen,
                                      ConfGenType const &end_conf_gen) const {
    return GTaskWorker<ConfGenType>(theta_, table_, m_, u_,
                                    start_conf_gen, end_conf_gen);
  }

 private:
  double theta_;
  Table *table_;
  uint32_t m_, u_;
};

#endif  // LDHELMET_PADE_G_TASK_H_
