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

#ifndef LDHELMET_TABLE_GEN_DEGREE_TASK_H_
#define LDHELMET_TABLE_GEN_DEGREE_TASK_H_

#include "common/conf_gen.h"
#include "common/vector_definitions.h"

// Performs work on given set of confs.
template<typename MarGenType, typename TaskType>
class DegreeTaskWorker {
 public:
  DegreeTaskWorker(TaskType const &task,
                     Vec8 *table,
                     double theta,
                     double rho,
                     MarGenType const &start_mar_gen,
                     MarGenType const &end_mar_gen)
      : task_(task),
        table_(table),
        theta_(theta),
        rho_(rho),
        start_mar_gen_(start_mar_gen),
        end_mar_gen_(end_mar_gen) { }

  void operator()() {
    for (MarGenType mar_gen = start_mar_gen_;
         mar_gen != end_mar_gen_;
         ++mar_gen) {
      assert(!mar_gen.end());
      Mar mar = *mar_gen;
      task_(table_, theta_, rho_,
           mar.a_mar_, mar.A_mar_, mar.b_mar_, mar.B_mar_);
    }
  }

 private:
  TaskType task_;
  Vec8 *table_;
  double theta_, rho_;
  MarGenType start_mar_gen_, end_mar_gen_;
};

// Contains parameters for work.
// Produces workers when given start and end mar iterators.
template<typename TaskType>
class DegreeTaskMaster {
 public:
  DegreeTaskMaster(TaskType const &task,
                     Vec8 *table,
                     double theta,
                     double rho)
      : task_(task),
        table_(table),
        theta_(theta),
        rho_(rho) {}

  template<typename MarGenType>
  DegreeTaskWorker<MarGenType, TaskType> operator()(
      MarGenType const &start_mar_gen,
      MarGenType const &end_mar_gen) const {
    return DegreeTaskWorker<MarGenType, TaskType>(
        task_, table_, theta_, rho_, start_mar_gen, end_mar_gen);
  }

 private:
  TaskType task_;
  Vec8 *table_;
  double theta_, rho_;
};

#endif  // LDHELMET_TABLE_GEN_DEGREE_TASK_H_
