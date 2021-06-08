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

#ifndef LDHELMET_TABLE_GEN_SOLVE_DEGREE_INL_H_
#define LDHELMET_TABLE_GEN_SOLVE_DEGREE_INL_H_

#include <stddef.h>
#include <stdint.h>

#include "common/conf.h"
#include "common/conf_gen.h"
#include "common/mar.h"
#include "common/predicates.h"
#include "common/threading.h"
#include "common/vector_definitions.h"
#include "table_gen/degree_task.h"
#include "table_gen/solve_scc.h"

template<class GenType>
size_t AssignBaseCasesMT(uint32_t num_threads,
                         uint64_t num_partitions,
                         double theta,
                         double rho,
                         GenType mar_gen,
                         Vec8 *table) {
  return RunThreaded(
    num_threads,
    DegreeTaskMaster<void(*)(Vec8 *, double, double, uint32_t,
                               uint32_t, uint32_t, uint32_t)>(
        &AssignBaseCases, table, theta, rho),
    PartitionTask(num_partitions, mar_gen));
}

template<class GenType>
size_t SolveRecursionMT(uint32_t num_threads,
                        uint64_t num_partitions,
                        double theta,
                        double rho,
                        GenType mar_gen,
                        Vec8 *table) {
  return RunThreaded(
    num_threads,
    DegreeTaskMaster<void(*)(Vec8 *, double, double, uint32_t,
                               uint32_t, uint32_t, uint32_t)>(
        &SolveSCC, table, theta, rho),
    PartitionTask(num_partitions, mar_gen));
}

template<class Pred>
void SolveDegree(double theta,
                 double rho,
                 Vec8 *table,
                 uint32_t degree,
                 uint32_t max_degree,
                 Pred const &pred,
                 size_t num_threads) {
  // A given degree can at most depend on itself, (degree - 1)
  // and (degree - 2).
  //
  // (degree-1) and (degree-2) should already be computed.
  // degree must be at least 2.

  // Check table allocation.
  if (table->size() < degree + 1) {
    fprintf(stderr,
            "Error: Table has not been allocated properly.\n");
    std::exit(1);
  }

  // Check that degree is at least 2.
  if (degree < 2) {
    fprintf(stderr, "Error: Degree must be at least 2.\n");
    std::exit(1);
  }

  // Partially validate that the table has the required information.
  // If degree == 2, then no prior information is necessary.
  // If degree is 3, then need one degree below.
  if (degree == 3 && (*table)[degree - 1].size() == 0) {
    fprintf(stderr,
            "Error: Degree 3 requires degree 2 in the table.\n");
  }

  // If degree > 3, then need at least two degrees below.
  if (degree > 3 && ((*table)[degree - 1].size() == 0 ||
                     (*table)[degree - 2].size() == 0)) {
    fprintf(stderr,
            "Error: Degree %d requires degrees %d and %d in the table.\n",
            static_cast<int>(degree),
            static_cast<int>(degree) - 1,
            static_cast<int>(degree) - 2);
    std::exit(1);
  }

  AllocateMemoryToTable(table, degree);

  // *** Begin solving SCCs ***
  // Solve in the following order:
  // 1. Base cases (i.e. configurations that have one ancestral allele
  //    remaining for one of the loci).
  // 2. Configurations that allow 0 mutations in current degree.
  // 3. Configurations that allow 1 mutation in current degree.
  // 4. Configurations that allow 2 mutations in current degree.

  uint64_t num_sccs = 0;
  uint64_t num_partitions = num_threads;

  // Base cases.
  {
    GenPred<MarGenBaseA, BaseCaseAPred<DecideToSolve<Pred> > > mar_gen(
      (MarGenBaseA(degree)),
      (BaseCaseAPred<DecideToSolve<Pred> >
        (DecideToSolve<Pred>(pred))));
    num_sccs += AssignBaseCasesMT(num_threads,
                                  num_partitions,
                                  theta,
                                  rho,
                                  mar_gen,
                                  table);
  }

  {
    GenPred<MarGenBaseB, BaseCaseBPred<DecideToSolve<Pred> > > mar_gen(
      (MarGenBaseB(degree)),
      (BaseCaseBPred<DecideToSolve<Pred> >
        (DecideToSolve<Pred>(pred))));
    num_sccs += AssignBaseCasesMT(num_threads,
                                  num_partitions,
                                  theta,
                                  rho,
                                  mar_gen,
                                  table);
  }

  // Inductive cases.
  // No mutation.
  {
    GenPred<MarGen, NoMutPred<DecideToSolve<Pred> > > mar_gen(
      (MarGen(degree)),
      (NoMutPred<DecideToSolve<Pred> >(DecideToSolve<Pred>(pred))));
    num_sccs += SolveRecursionMT(num_threads,
                                 num_partitions,
                                 theta,
                                 rho,
                                 mar_gen,
                                 table);
  }

  // Mutation at locus 1.
  {
    GenPred<MarGenMutA, MutAPred<DecideToSolve<Pred> > > mar_gen(
      (MarGenMutA(degree)),
      (MutAPred<DecideToSolve<Pred> >(DecideToSolve<Pred>(pred))));
    num_sccs += SolveRecursionMT(num_threads,
                                 num_partitions,
                                 theta,
                                 rho,
                                 mar_gen,
                                 table);
  }

  // Mutation at locus 2.
  {
    GenPred<MarGenMutB, MutBPred<DecideToSolve<Pred> > > mar_gen(
      (MarGenMutB(degree)),
      (MutBPred<DecideToSolve<Pred> >(DecideToSolve<Pred>(pred))));
    num_sccs += SolveRecursionMT(num_threads,
                                 num_partitions,
                                 theta,
                                 rho,
                                 mar_gen,
                                 table);
  }

  // Mutations at loci 1 and 2.
  {
    GenPred<MarGenMutAB, MutABPred<DecideToSolve<Pred> > > mar_gen(
      (MarGenMutAB(degree)),
      (MutABPred<DecideToSolve<Pred> >(DecideToSolve<Pred>(pred))));
    num_sccs += SolveRecursionMT(num_threads,
                                 num_partitions,
                                 theta,
                                 rho,
                                 mar_gen,
                                 table);
  }

  // Free memory in table for degree - 2, which is no longer needed by
  // subsequent degrees.
  if (degree >= 4) {
    Vec7().swap((*table)[degree - 2]);
  }
}

#endif  // LDHELMET_TABLE_GEN_SOLVE_DEGREE_INL_H_
