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

#include "common/conf.h"

#include <stdint.h>

#include <string>

#include <boost/lexical_cast.hpp>

Conf::Conf() : a_(0), A_(0), b_(0), B_(0), ab_(0), aB_(0), Ab_(0), AB_(0) { }

Conf::Conf(uint8_t a, uint8_t A, uint8_t b, uint8_t B,
           uint8_t ab, uint8_t aB, uint8_t Ab, uint8_t AB)
    : a_(a), A_(A), b_(b), B_(B), ab_(ab), aB_(aB), Ab_(Ab), AB_(AB) { }

Conf::Conf(BinaryRep const &binary_rep) {
  a_  = (binary_rep.rep() >> 0)  & 0xffL;
  A_  = (binary_rep.rep() >> 8)  & 0xffL;
  b_  = (binary_rep.rep() >> 16) & 0xffL;
  B_  = (binary_rep.rep() >> 24) & 0xffL;
  ab_ = (binary_rep.rep() >> 32) & 0xffL;
  aB_ = (binary_rep.rep() >> 40) & 0xffL;
  Ab_ = (binary_rep.rep() >> 48) & 0xffL;
  AB_ = (binary_rep.rep() >> 56) & 0xffL;
}

std::string Conf::GetString() const {
  return boost::lexical_cast<std::string>(static_cast<uint32_t>(a_))  + " " +
         boost::lexical_cast<std::string>(static_cast<uint32_t>(A_))  + " " +
         boost::lexical_cast<std::string>(static_cast<uint32_t>(b_))  + " " +
         boost::lexical_cast<std::string>(static_cast<uint32_t>(B_))  + " " +
         boost::lexical_cast<std::string>(static_cast<uint32_t>(ab_)) + " " +
         boost::lexical_cast<std::string>(static_cast<uint32_t>(aB_)) + " " +
         boost::lexical_cast<std::string>(static_cast<uint32_t>(Ab_)) + " " +
         boost::lexical_cast<std::string>(static_cast<uint32_t>(AB_));
}

Conf Conf::Adda(int i, int quantity) const {
  Conf tmp_conf = *this;
  if (i == 0) {
    tmp_conf.a_ += quantity;
  } else {
    tmp_conf.A_ += quantity;
  }
  return tmp_conf;
}

Conf Conf::Addb(int j, int quantity) const {
  Conf tmp_conf = *this;
  if (j == 0) {
    tmp_conf.b_ += quantity;
  } else {
    tmp_conf.B_ += quantity;
  }
  return tmp_conf;
}

Conf Conf::Addr(int i, int j, int quantity) const {
  Conf tmp_conf = *this;
  if (i == 0) {
    if (j == 0) {
      tmp_conf.ab_ += quantity;
    } else {
      tmp_conf.aB_ += quantity;
    }
  } else {
    if (j == 0) {
      tmp_conf.Ab_ += quantity;
    } else {
      tmp_conf.AB_ += quantity;
    }
  }
  return tmp_conf;
}

bool Conf::operator<(Conf const &conf2) const {
  Conf const &conf1 = *this;
  uint32_t level1 = conf1.a_ + conf1.A_ + conf1.b_ + conf1.B_
              + 2*(conf1.ab_ + conf1.aB_ + conf1.Ab_ + conf1.AB_);
  uint32_t level2 = conf2.a_ + conf2.A_ + conf2.b_ + conf2.B_
              + 2*(conf2.ab_ + conf2.aB_ + conf2.Ab_ + conf2.AB_);
  if (level1 != level2) {
    return level1 < level2;
  } else if (conf1.a_  != conf2.a_) {
    return conf1.a_  <  conf2.a_;
  } else if (conf1.A_  != conf2.A_) {
    return conf1.A_  <  conf2.A_;
  } else if (conf1.b_  != conf2.b_) {
    return conf1.b_  <  conf2.b_;
  } else if (conf1.B_  != conf2.B_) {
    return conf1.B_  <  conf2.B_;
  } else if (conf1.ab_ != conf2.ab_) {
    return conf1.ab_ <  conf2.ab_;
  } else if (conf1.aB_ != conf2.aB_) {
    return conf1.aB_ <  conf2.aB_;
  } else if (conf1.Ab_ != conf2.Ab_) {
    return conf1.Ab_ <  conf2.Ab_;
  } else {
    return conf1.AB_ < conf2.AB_;
  }
}

Conf::BinaryRep::BinaryRep() : rep_(0) { }

Conf::BinaryRep::BinaryRep(Conf const &conf) : rep_(0) {
  rep_ += (static_cast<uint64_t>(conf.a_))  << 0;
  rep_ += (static_cast<uint64_t>(conf.A_))  << 8;
  rep_ += (static_cast<uint64_t>(conf.b_))  << 16;
  rep_ += (static_cast<uint64_t>(conf.B_))  << 24;
  rep_ += (static_cast<uint64_t>(conf.ab_)) << 32;
  rep_ += (static_cast<uint64_t>(conf.aB_)) << 40;
  rep_ += (static_cast<uint64_t>(conf.Ab_)) << 48;
  rep_ += (static_cast<uint64_t>(conf.AB_)) << 56;
}
