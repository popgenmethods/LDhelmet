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

#ifndef LDHELMET_PADE_Q_TASK_H_
#define LDHELMET_PADE_Q_TASK_H_

#include <stdint.h>

#include <vector>

#include "common/conf_gen.h"
#include "pade/coeff.h"

// Performs work on given set of confs.
template<typename ConfIDGenType>
class QTaskWorker {
 public:
  QTaskWorker(Table *table,
              std::vector<Conf> const &conf_list,
              std::vector<std::vector<double> > &q_values,
              uint32_t cur_level,
              uint32_t m,
              ConfIDGenType const &start_conf_id_gen,
              ConfIDGenType const &end_conf_id_gen)
      : table_(table),
        conf_list_(conf_list),
        q_values_(q_values),
        cur_level_(cur_level),
        m_(m),
        start_conf_id_gen_(start_conf_id_gen),
        end_conf_id_gen_(end_conf_id_gen) { }

  void operator()() {
    for (ConfIDGenType conf_id_gen = start_conf_id_gen_;
         conf_id_gen != end_conf_id_gen_;
         ++conf_id_gen) {
      assert(!conf_id_gen.end());
      q_values_[*conf_id_gen][cur_level_/2] +=
          ComputeQValueIncrementally(*table_,
                                     cur_level_,
                                     m_,
                                     conf_list_[*conf_id_gen]);
    }
  }

 private:
  Table *table_;
  std::vector<Conf> const &conf_list_;
  std::vector<std::vector<double> > &q_values_;
  uint32_t cur_level_;
  uint32_t m_;
  ConfIDGenType start_conf_id_gen_, end_conf_id_gen_;
};

// Contains parameters for work.
// Produces workers when given start and end conf iterators.
class QTaskMaster {
 public:
  QTaskMaster(Table *table,
              std::vector<Conf> const &conf_list,
              std::vector<std::vector<double> > &q_values,
              uint32_t cur_level,
              uint32_t m)
      : table_(table),
        conf_list_(conf_list),
        q_values_(q_values),
        cur_level_(cur_level),
        m_(m) { }

  template<typename ConfIDGenType>
  QTaskWorker<ConfIDGenType> operator()(
      ConfIDGenType const &start_conf_id_gen,
      ConfIDGenType const &end_conf_id_gen) const {
    return QTaskWorker<ConfIDGenType>(table_,
                                      conf_list_,
                                      q_values_,
                                      cur_level_,
                                      m_,
                                      start_conf_id_gen,
                                      end_conf_id_gen);
  }

 private:
  Table *table_;
  std::vector<Conf> const &conf_list_;
  std::vector<std::vector<double> > &q_values_;
  uint32_t cur_level_;
  uint32_t m_;
};

#endif  // LDHELMET_PADE_Q_TASK_H_
