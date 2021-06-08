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

#ifndef LDHELMET_RJMCMC_RAN_NUM_GEN_H_
#define LDHELMET_RJMCMC_RAN_NUM_GEN_H_

// Header file implementing a random number generator class.

#include <stdio.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

// Default seed for random number generator.
static boost::mt19937::result_type const kDefaultSeed = 5489;

class RanNumGen {
 public:
  typedef boost::mt19937 RngType;
  typedef RngType::result_type SeedType;

  typedef boost::variate_generator<RngType&, boost::uniform_real<> >
    UniformGenBaseType;
  typedef boost::variate_generator<RngType&, boost::uniform_int<> >
    UniformIntGenType;

  typedef RngType RandomNumberGenType;

  // Continuous uniform distribution from lower_bound to upper_bound,
  // exclusive.
  class UniformGenType {
   public:
    UniformGenType(RanNumGen *rng,
                   double lower_bound,
                   double upper_bound)
        : uniform_(rng->GetUniformGenBase(lower_bound, upper_bound)),
          lower_bound_(lower_bound),
          upper_bound_(upper_bound) {
      if (lower_bound_ >= upper_bound_) {
        fprintf(stderr,
                "Error: UniformGenType: lower_bound >= upper_bound.\n");
        std::exit(1);
      }
    }

    double operator()() {
      // Rejection sampling to exclude boundaries.
      double u = uniform_();
      while (u == lower_bound_ || u == upper_bound_) {
        u = uniform_();
      }
      return u;
    }

    private:
      UniformGenBaseType uniform_;
      double lower_bound_, upper_bound_;
  };

  RanNumGen() : seed_(kDefaultSeed), gen_(seed_) { }
  explicit RanNumGen(SeedType seed) : seed_(seed), gen_(seed_) { }

  // Continuous uniform distribution from lower_bound to upper_bound,
  // inclusive.
  UniformGenBaseType GetUniformGenBase(double lower_bound,
                                       double upper_bound) {
    boost::uniform_real<> uniform_dist(lower_bound, upper_bound);
    UniformGenBaseType uniform(gen_, uniform_dist);
    return uniform;
  }

  // Discrete uniform distribution from lower_bound to upper_bound,
  // inclusive.
  UniformIntGenType GetUniformIntGen(int lower_bound,
                                     int upper_bound) {
    boost::uniform_int<> uniform_int_dist(lower_bound, upper_bound);
    boost::variate_generator<RngType&, boost::uniform_int<> >
      uniform_int(gen_, uniform_int_dist);
    return uniform_int;
  }

  RandomNumberGenType& GetRandomNumberGen() {
    return gen_;
  }

  UniformGenType GetUniformGen(double lower_bound,
                               double upper_bound) {
    return UniformGenType(this, lower_bound, upper_bound);
  }

 private:
  SeedType seed_;
  RngType gen_;
};

#endif  // LDHELMET_RJMCMC_RAN_NUM_GEN_H_
