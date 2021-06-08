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

#include "common/conf_gen.h"

#include <stdint.h>

#include <algorithm>
#include <cassert>

#include "common/conf.h"
#include "common/mar.h"
#include "common/predicates.h"

MarGenBase::MarGenBase(uint8_t degree) : end_(false), degree_(degree) { }

MarGenBaseA::MarGenBaseA() : MarGenBase(0) { }

MarGenBaseA::MarGenBaseA(uint8_t degree) : MarGenBase(degree) {
  if (degree_ == 0) {
    return;
  }
  mar_.a_mar_ = 1;
  max_a_mar_  = 1;

  mar_.A_mar_ = 0;
  max_A_mar_  = 0;

  mar_.b_mar_ = 0;
  max_b_mar_  = degree_ - 1;

  mar_.B_mar_ = degree_ - 1 - mar_.b_mar_;
  max_B_mar_  = degree_ - 1 - mar_.b_mar_;

  assert(mar_.ComputeDegree() == degree_);
}

MarGenBaseA &MarGenBaseA::operator++() {
  if (end()) {
    return *this;
  }
  ++mar_.b_mar_;
  if (mar_.b_mar_ > max_b_mar_) {
    end_ = true;
    return *this;
  }

  mar_.B_mar_ = degree_ - 1 - mar_.b_mar_;

  return *this;
}

MarGenBaseB::MarGenBaseB() : MarGenBase(0) { }

MarGenBaseB::MarGenBaseB(uint8_t degree) : MarGenBase(degree) {
  if (degree_ == 0) {
    return;
  }

  mar_.a_mar_ = 0;
  max_a_mar_  = degree_ - 1;

  mar_.A_mar_ = degree_ - 1 - mar_.a_mar_;
  max_A_mar_  = degree_ - 1 - mar_.a_mar_;

  mar_.b_mar_ = 1;
  max_b_mar_  = 1;

  mar_.B_mar_ = 0;
  max_B_mar_  = 0;

  assert(mar_.ComputeDegree() == degree_);
}

MarGenBaseB &MarGenBaseB::operator++() {
  if (end()) {
    return *this;
  }
  ++mar_.a_mar_;
  if (mar_.a_mar_ > max_a_mar_) {
      end_ = true;
      return *this;
  }
  mar_.A_mar_ = degree_ - 1 - mar_.a_mar_;

  return *this;
}

MarGenMutA::MarGenMutA() : MarGenBase(0) { }

MarGenMutA::MarGenMutA(uint8_t degree) : MarGenBase(degree) {
  if (degree_ == 0) {
    return;
  }

  assert(degree > 0);

  mar_.a_mar_ = 0;
  max_a_mar_  = degree_ - 1;

  mar_.A_mar_ = 1;
  max_A_mar_  = 1;

  mar_.b_mar_ = 0;
  max_b_mar_  = degree_ - mar_.a_mar_ - mar_.A_mar_;

  mar_.B_mar_ = degree_ - mar_.a_mar_ - mar_.A_mar_ - mar_.b_mar_;
  max_B_mar_  = degree_ - mar_.a_mar_ - mar_.A_mar_ - mar_.b_mar_;

  assert(mar_.ComputeDegree() == degree_);
}

MarGenMutA &MarGenMutA::operator++() {
    if (end()) {
      return *this;
    }
    ++mar_.b_mar_;
    if (mar_.b_mar_ > max_b_mar_) {
      ++mar_.a_mar_;
      if (mar_.a_mar_ > max_a_mar_) {
          end_ = true;
          return *this;
      }
      max_b_mar_ = degree_ - mar_.a_mar_ - mar_.A_mar_;
      mar_.b_mar_ = 0;
    }
    mar_.B_mar_ = degree_ - mar_.a_mar_ - mar_.A_mar_ - mar_.b_mar_;

    return *this;
}

MarGenMutB::MarGenMutB() : MarGenBase(0) { }

MarGenMutB::MarGenMutB(uint8_t degree) : MarGenBase(degree) {
  if (degree == 0) {
    return;
  }
  mar_.a_mar_ = 0;
  max_a_mar_  = degree_ - 1;

  mar_.A_mar_ = 0;
  max_A_mar_  = degree_ - mar_.a_mar_ - 1;

  mar_.b_mar_ = degree_ - mar_.a_mar_ - mar_.A_mar_ - 1;
  max_b_mar_  = degree_ - mar_.a_mar_ - mar_.A_mar_ - 1;

  mar_.B_mar_ = 1;
  max_B_mar_  = 1;

  assert(mar_.ComputeDegree() == degree_);
}

MarGenMutB &MarGenMutB::operator++() {
  if (end()) {
    return *this;
  }
  ++mar_.A_mar_;
  if (mar_.A_mar_ > max_A_mar_) {
    ++mar_.a_mar_;
    if (mar_.a_mar_ > max_a_mar_) {
      end_ = true;
      return *this;
    }
    max_A_mar_ = degree_ - mar_.a_mar_ -1;
    mar_.A_mar_ = 0;
  }
  mar_.b_mar_ = degree_ - mar_.a_mar_ - mar_.A_mar_ - 1;

  return *this;
}

MarGenMutAB::MarGenMutAB() : MarGenBase(0) { }

MarGenMutAB::MarGenMutAB(uint8_t degree) : MarGenBase(degree) {
  if (degree_ == 0) {
    return;
  }
  mar_.a_mar_ = 0;
  max_a_mar_  = degree_ - 2;

  mar_.A_mar_ = 1;
  max_A_mar_  = 1;

  mar_.b_mar_ = degree_ - mar_.a_mar_ - 2;
  max_b_mar_  = degree_ - mar_.a_mar_ - 2;

  mar_.B_mar_ = 1;
  max_B_mar_  = 1;

  assert(mar_.ComputeDegree() == degree_);
}

MarGenMutAB &MarGenMutAB::operator++() {
  if (end()) {
    return *this;
  }
  ++mar_.a_mar_;
  if (mar_.a_mar_ > max_a_mar_) {
      end_ = true;
      return *this;
  }
  mar_.b_mar_ = degree_ - mar_.a_mar_ - 2;

  return *this;
}

MarGen::MarGen() : MarGenBase(0) { }

MarGen::MarGen(uint8_t degree) : MarGenBase(degree) {
  if (degree_ == 0) {
    return;
  }

  mar_.a_mar_ = 0;
  mar_.A_mar_ = 0;
  mar_.b_mar_ = 0;
  mar_.B_mar_ = degree_ - mar_.a_mar_ - mar_.A_mar_ - mar_.b_mar_;

  max_a_mar_ = degree_;
  max_A_mar_ = max_a_mar_ - mar_.a_mar_;
  max_b_mar_ = max_A_mar_ - mar_.A_mar_;

  assert(mar_.ComputeDegree() == degree_);
}

MarGen &MarGen::operator++() {
  if (end()) {
    return *this;
  }
  ++mar_.b_mar_;
  if (mar_.b_mar_ > max_b_mar_) {
    ++mar_.A_mar_;
    if (mar_.A_mar_ > max_A_mar_) {
      ++mar_.a_mar_;
      if (mar_.a_mar_ > max_a_mar_) {
        end_ = true;
        return *this;
      }
      max_A_mar_ = max_a_mar_ - mar_.a_mar_;
      mar_.A_mar_ = 0;
    }
    max_b_mar_ = max_A_mar_ - mar_.A_mar_;
    mar_.b_mar_ = 0;
  }
  mar_.B_mar_ = degree_ - mar_.a_mar_ - mar_.A_mar_ - mar_.b_mar_;

  return *this;
}

ConfGenM::ConfGenM()
    : end_(false), max_degree_(0), m_(0) {
  Init(max_degree_, m_);
}

ConfGenM::ConfGenM(uint8_t max_degree, uint8_t m)
    : end_(false), max_degree_(max_degree), m_(m) {
  Init(max_degree_, m_);
}

void ConfGenM::Init(uint8_t max_degree, uint8_t m) {
  assert(max_degree_ >= 2*m);

  degree_ = std::max(1, 2*m);

  conf_.ab_ = 0;

  conf_.aB_ = 0;
  max_aB_   = m_ - conf_.ab_;

  conf_.Ab_ = 0;
  max_Ab_   = m_ - conf_.ab_ - conf_.aB_;

  conf_.AB_ = m_ - conf_.ab_ - conf_.aB_ - conf_.Ab_;

  conf_.a_  = 0;
  max_a_    = degree_ - 2*m_;

  conf_.A_  = 0;
  max_A_    = max_a_ - conf_.a_;

  conf_.b_  = 0;
  max_b_    = max_a_ - conf_.a_ - conf_.A_;

  conf_.B_  = max_a_ - conf_.a_ - conf_.A_ - conf_.b_;

  assert(conf_.ab_ + conf_.aB_ + conf_.Ab_ + conf_.AB_ == m_);
  assert(conf_.a_ + conf_.A_ + conf_.b_ + conf_.B_
       + 2*(conf_.ab_ + conf_.aB_ + conf_.Ab_ + conf_.AB_) == degree_);
}

ConfGenM &ConfGenM::operator++() {
  assert(!end());
  assert(conf_.ab_ + conf_.aB_ + conf_.Ab_ + conf_.AB_ == m_);
  assert(conf_.a_ + conf_.A_ + conf_.b_ + conf_.B_
       + 2*(conf_.ab_ + conf_.aB_ + conf_.Ab_ + conf_.AB_) == degree_);

  if (end()) {
    return *this;
  }

  ++conf_.b_;
  if (conf_.b_ > max_b_) {
    ++conf_.A_;
    if (conf_.A_ > max_A_) {
      ++conf_.a_;
      if (conf_.a_ > max_a_) {
        ++conf_.Ab_;
        if (conf_.Ab_ > max_Ab_) {
          ++conf_.aB_;
          if (conf_.aB_ > max_aB_) {
            ++conf_.ab_;
            if (conf_.ab_ > m_) {
              ++degree_;
              if (degree_ > max_degree_) {
                  end_ = true;
                  return *this;
              }
              conf_.ab_ = 0;
            }
            max_aB_ = m_ - conf_.ab_;
            conf_.aB_ = 0;
          }
          max_Ab_ = m_ - conf_.ab_ - conf_.aB_;
          conf_.Ab_ = 0;
        }
        conf_.AB_ = m_ - conf_.ab_ - conf_.aB_ - conf_.Ab_;
        max_a_ = degree_ - 2*m_;
        conf_.a_ = 0;
      }
      max_A_ = max_a_ - conf_.a_;
      conf_.A_ = 0;
    }
    max_b_ = max_A_ - conf_.A_;
    conf_.b_ = 0;
  }
  conf_.B_ = max_b_ - conf_.b_;

  return *this;
}

ConfGenMC::ConfGenMC(Conf const &c, uint8_t m)
    : end_(false), c_(c), m_(m) {
  degree_ = c_.ComputeDegree();
  assert(degree_ >= 2*m && degree_ >= 1);

  uint8_t cTot = c_.ab_ + c_.aB_ + c_.Ab_ + c_.AB_;
  assert(m_ <= cTot);

  max_ab_ = std::min(m_, c_.ab_);
  conf_.ab_ = std::max(0, m_ - c_.aB_ - c_.Ab_ - c_.AB_);

  max_aB_ = std::min(static_cast<uint8_t>(m_ - conf_.ab_), c_.aB_);
  conf_.aB_ = std::max(0, (m_ - conf_.ab_) - c_.Ab_ - c_.AB_);

  max_Ab_ = std::min(static_cast<uint8_t>(m_ - conf_.ab_ - conf_.aB_), c_.Ab_);
  conf_.Ab_ = std::max(0, (m_ - conf_.ab_ - conf_.aB_) - c_.AB_);

  conf_.AB_ = m_ - conf_.ab_ - conf_.aB_ - conf_.Ab_;
  assert(conf_.AB_ <= c_.AB_);

  conf_.a_ = c_.a_ + (c_.ab_ + c_.aB_) - (conf_.ab_ + conf_.aB_);
  conf_.A_ = c_.A_ + (c_.Ab_ + c_.AB_) - (conf_.Ab_ + conf_.AB_);
  conf_.b_ = c_.b_ + (c_.ab_ + c_.Ab_) - (conf_.ab_ + conf_.Ab_);
  conf_.B_ = c_.B_ + (c_.aB_ + c_.AB_) - (conf_.aB_ + conf_.AB_);
}

ConfGenMC &ConfGenMC::operator++() {
  assert(!end());
  assert(conf_.ab_ + conf_.aB_ + conf_.Ab_ + conf_.AB_ == m_);
  assert(conf_.ComputeDegree() == degree_);

  if (end()) {
    return *this;
  }

  ++conf_.Ab_;
  if (conf_.Ab_ > max_Ab_) {
    ++conf_.aB_;
    if (conf_.aB_ > max_aB_) {
      ++conf_.ab_;
      if (conf_.ab_ > max_ab_) {
        end_ = true;
        return *this;
      }
      assert(end() || (conf_.ab_ <= c_.ab_ && conf_.ab_ >= 0));

      max_aB_ = std::min(static_cast<uint8_t>(m_ - conf_.ab_), c_.aB_);
      conf_.aB_ = std::max(0, (m_ - conf_.ab_) - c_.Ab_ - c_.AB_);

      assert(end() || (conf_.aB_ <= c_.aB_ && conf_.aB_ >= 0));
    }
    max_Ab_ = std::min(static_cast<uint8_t>(m_ - conf_.ab_ - conf_.aB_),
                       c_.Ab_);
    conf_.Ab_ = std::max(0, (m_ - conf_.ab_ - conf_.aB_) - c_.AB_);
    assert(end() || (conf_.Ab_ <= c_.Ab_ && conf_.Ab_ >= 0));
  }
  conf_.AB_ = m_ - conf_.ab_ - conf_.aB_ - conf_.Ab_;
  assert(end() || (conf_.AB_ <= c_.AB_ && conf_.AB_ >= 0));

  conf_.a_ = c_.a_ + (c_.ab_ + c_.aB_) - (conf_.ab_ + conf_.aB_);
  conf_.A_ = c_.A_ + (c_.Ab_ + c_.AB_) - (conf_.Ab_ + conf_.AB_);
  conf_.b_ = c_.b_ + (c_.ab_ + c_.Ab_) - (conf_.ab_ + conf_.Ab_);
  conf_.B_ = c_.B_ + (c_.aB_ + c_.AB_) - (conf_.aB_ + conf_.AB_);
  assert(end() || conf_.a_ >= 0);
  assert(end() || conf_.A_ >= 0);
  assert(end() || conf_.b_ >= 0);
  assert(end() || conf_.B_ >= 0);

  return *this;
}

ConfIDGen::ConfIDGen()
    : end_(false),
      start_id_(0),
      end_id_(0),
      cur_conf_id_(0) {
  if (cur_conf_id_ >= end_id_) {
    end_ = true;
  }
}

ConfIDGen::ConfIDGen(uint64_t start_id, uint64_t end_id)
    : end_(false),
      start_id_(start_id),
      end_id_(end_id),
      cur_conf_id_(start_id) {
  if (cur_conf_id_ >= end_id_) {
    end_ = true;
  }
}

ConfIDGen &ConfIDGen::operator++() {
  if (end()) {
    return *this;
  } else {
    ++cur_conf_id_;
    if (cur_conf_id_ >= end_id_) {
      end_ = true;
    }
    return *this;
  }
}
