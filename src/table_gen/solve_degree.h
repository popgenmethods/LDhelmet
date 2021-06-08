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

#ifndef LDHELMET_TABLE_GEN_SOLVE_DEGREE_H_
#define LDHELMET_TABLE_GEN_SOLVE_DEGREE_H_

#include <stddef.h>
#include <stdint.h>

#include "common/conf.h"
#include "common/mar.h"
#include "common/predicates.h"
#include "common/vector_definitions.h"

// Assign base cases using the single locus formula.
void AssignBaseCases(Vec8 *table,
                     double theta,
                     double rho,
                     uint32_t a_mar,
                     uint32_t A_mar,
                     uint32_t b_mar,
                     uint32_t B_mar);

template<class GenType>
size_t AssignBaseCasesMT(uint32_t num_threads,
                        uint64_t num_partitions,
                        double theta,
                        double rho,
                        GenType mar_gen,
                        Vec8* table);

template<class GenType>
size_t SolveRecursionMT(uint32_t num_threads,
                    uint64_t num_partitions,
                    double theta,
                    double rho,
                    GenType mar_gen,
                    Vec8* table);

// Given a degree, solveDegree() will find the probabilities of every
// configuration in that degree.
// A configuration is in a given degree if the number of defined alleles is
// equal to the degree.
// i.e. the degree of a configuration is equal to a+b+2*c.
// table contains the probabilities for some number of degrees (usually 3
// degrees at once)
template<class Pred>
void SolveDegree(double theta,
                 double rho,
                 Vec8 *table,
                 uint32_t degree,
                 uint32_t max_degree,
                 Pred const &pred,
                 size_t num_threads);

#endif  // LDHELMET_TABLE_GEN_SOLVE_DEGREE_H_
