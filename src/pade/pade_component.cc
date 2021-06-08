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

#include <stdint.h>
#include <stdio.h>

#include <vector>
#include <cstdlib>

#include <boost/tuple/tuple.hpp>

#include "common/conf_gen.h"
#include "common/read_confs.h"
#include "common/threading.h"
#include "common/version_number.h"
#include "pade/coeff.h"
#include "pade/compute_g.h"
#include "pade/g_task.h"
#include "pade/memory_manager.h"
#include "pade/memory_manager-inl.h"
#include "pade/pade_options.h"
#include "pade/q_task.h"
#include "pade/subtable.h"

int PadeComponent(std::string const base_command, int argc, char **argv) {
  uint64_t version_number = PADE_VERSION_OUTPUT;
  uint64_t version_bit_string = PADE_SALT + version_number;

  std::string version_string = PADE_VERSION_STRING;
  CmdLineOptionsPade cmd_line_options(base_command, argc, argv, version_string);

  if (!cmd_line_options.success()) {
    exit(1);
  }

  printf("%s\n", ShowCommandLineOptions(argc, argv).c_str());

  uint32_t const &num_threads = cmd_line_options.num_threads_;
  uint64_t const &num_coeffs = cmd_line_options.num_coeffs_;

  printf("Number of threads: %d.\n",
         static_cast<int>(cmd_line_options.num_threads_));
  assert(num_threads > 0);

  printf("Number of coefficients: %d.\n", static_cast<int>(num_coeffs));

  std::string const &conf_file = cmd_line_options.conf_file_;
  std::string const &output_file = cmd_line_options.output_file_;

  double const &theta = cmd_line_options.theta_;
  double const &defect_threshold = cmd_line_options.defect_threshold_;

  if (num_threads < 1) {
    fprintf(stderr,
            "Number of threads must be at least 1.\n");
    std::exit(1);
  }

  if (num_coeffs < 2) {
    fprintf(stderr,
            "Number of coefficients must be at least 2.\n");
    std::exit(1);
  }

  std::vector<Conf> conf_list = LoadConfigurations(conf_file);

  uint32_t max_degree, max_locus;
  boost::tie(boost::tuples::ignore,
             max_degree,
             boost::tuples::ignore,
             max_locus) = PreprocessConfs(conf_list);

  std::vector<std::vector<double> >
    q_values(conf_list.size(), std::vector<double>(num_coeffs, 0.0));

  // Global table: All threads share.
  assert(num_coeffs >= 2);
  Table table(2*num_coeffs - 1, 2*num_coeffs - 1);

  MemoryManager<GenPred<ConfGenM, SCCMutation> > memory_manager(&table,
                                                                max_degree);

  printf("Computing Pade coefficients.\n");

  uint64_t num_confs = 0;
  for (uint32_t cur_level = 0; cur_level < 2*num_coeffs; cur_level += 2) {
    printf("Level: %d\n", static_cast<int>(cur_level));

    // Note: max_degree/2 >= c_total, but c_total might be a tighter bound.

    for (uint32_t m = 0; m <= std::min(cur_level, max_degree / 2); ++m) {
      uint32_t const u = cur_level - m;  // Skip m + u = l where l is odd.

      printf("(m:%d, u:%d)", static_cast<int>(m), static_cast<int>(u));
      fflush(stdout);

      // Compute g_u^(m).
      {
        GenPred<ConfGenM, SCCMutation> const conf_gen(ConfGenM(max_degree, m),
                                                     SCCMutation(max_locus));

        // Read in necessary parts of previous antidiagonal.
        if (m == 0) {
          memory_manager.LoadAllElementsNeeded(m, u, conf_gen);
        } else {
          memory_manager.LoadElementsIncrementally(m, u, conf_gen);
        }

        // Allocate space to table for current m and u.
        table.AllocateMemory(max_degree, m, u);

        size_t num_confs_gen =
          RunThreaded(num_threads,
                      GTaskMaster(theta, &table, m, u),
                      PartitionTask(num_threads, conf_gen));

        num_confs += num_confs_gen;

        // Compute q values.
        ConfIDGen const conf_id_gen(0, conf_list.size());
        RunThreaded(num_threads,
                    QTaskMaster(&table,
                                conf_list,
                                q_values,
                                cur_level,
                                m),
                    PartitionTask(num_threads, conf_id_gen));

        // Write g values to file for current m and u so that it can be
        // used for the next antidiagonal.

        memory_manager.WriteElementsToFile(m, u, conf_gen);

        // Free memory of current m and u.
        memory_manager.FreeMU(m, u);

        // Free memory on previous antidiagonal not needed for next m
        // and u.
        memory_manager.FreeElementsIncrementally(m, u);
      }
    }

    printf("\n");
  }

  printf("Writing results to %s.\n", output_file.c_str());

  {
    // Write pade coefficients to file.
    FILE *fp = fopen(output_file.c_str(), "wb");
    if (fp == NULL) {
      fprintf(stderr,
              "Error opening output file for writing.\n");
      std::exit(1);
    }

    // Write header information.

    // Write version number.
    int num_written;
    num_written = fwrite(reinterpret_cast<char const *>(&version_bit_string),
                  sizeof(uint64_t), 1, fp);
    assert(num_written = 1);

    uint64_t num_conf_list = conf_list.size();
    num_written = fwrite(reinterpret_cast<char const *>(&num_conf_list),
                         sizeof(num_conf_list), 1, fp);
    assert(num_written = 1);
    num_written = fwrite(reinterpret_cast<char const *>(&theta),
                         sizeof(theta), 1, fp);
    assert(num_written = 1);
    num_written = fwrite(reinterpret_cast<char const *>(&num_coeffs),
                         sizeof(num_coeffs), 1, fp);
    assert(num_written = 1);
    num_written = fwrite(reinterpret_cast<char const *>(&defect_threshold),
                         sizeof(defect_threshold), 1, fp);
    assert(num_written = 1);

    // Write confs.
    for (std::vector<Conf>::const_iterator citer = conf_list.begin();
         citer != conf_list.end();
         ++citer) {
      Conf::BinaryRep::RepType rep = Conf::BinaryRep(*citer).rep();
      num_written = fwrite(reinterpret_cast<char const *>(&rep),
                           sizeof(rep), 1, fp);
      assert(num_written = 1);
    }

    // Write pade coefficients and roots to file.
    for (size_t conf_id = 0; conf_id < conf_list.size(); ++conf_id) {
      Coeffs coeffs = ComputeCoeffsForConf(num_coeffs,
                                           q_values[conf_id],
                                           conf_list[conf_id]);
      for (size_t coeff_id = 0; coeff_id < num_coeffs; ++coeff_id) {
        double pade_coeff = coeffs.cont_frac[coeff_id];
        num_written = fwrite(reinterpret_cast<char const *>(&pade_coeff),
                             sizeof(pade_coeff), 1, fp);
        assert(num_written = 1);
      }

      PadeRoots roots = ComputeRelevantRootsForConf(num_coeffs,
                                                    defect_threshold,
                                                    coeffs);

      // Combine roots from numerator and denominator.
      std::vector<std::vector<double> > conf_roots;
      std::vector<std::vector<double> > const &numerator = roots.numerator;
      std::vector<std::vector<double> > const &denominator = roots.denominator;
      size_t max_num_coeffs_with_roots = std::max(numerator.size(),
                                                  denominator.size());

      for (size_t coeff_root_id = 0;
           coeff_root_id < max_num_coeffs_with_roots;
           ++coeff_root_id) {
        conf_roots.push_back(std::vector<double>());
        if (numerator.size() > coeff_root_id) {
            for (size_t i = 0; i < numerator[coeff_root_id].size(); ++i) {
              conf_roots.back().push_back(numerator[coeff_root_id][i]);
              assert(numerator[coeff_root_id][i] >= 0);
            }
        }
        if (denominator.size() > coeff_root_id) {
          for (size_t i = 0; i < denominator[coeff_root_id].size(); ++i) {
            conf_roots.back().push_back(denominator[coeff_root_id][i]);
            assert(denominator[coeff_root_id][i] >= 0);
          }
        }
        assert(conf_roots.back().size() > 0);
      }

      // Write number of sets of coefficients with roots for current conf.
      assert(conf_roots.size() < 256);
      uint8_t num_coeffs_with_roots = static_cast<uint8_t>(conf_roots.size());
      num_written = fwrite(
          reinterpret_cast<char const *>(&num_coeffs_with_roots),
          sizeof(num_coeffs_with_roots), 1, fp);
      assert(num_written = 1);

      for (size_t coeff_root_id = 0;
         coeff_root_id < conf_roots.size();
         ++coeff_root_id) {
        // Write number of roots for given set of coefficients.
        assert(conf_roots[coeff_root_id].size() < 256);
        uint8_t num_roots_for_coeff =
          static_cast<uint8_t>(conf_roots[coeff_root_id].size());
        num_written = fwrite(
            reinterpret_cast<char const *>(&num_roots_for_coeff),
            sizeof(num_roots_for_coeff), 1, fp);
        assert(num_written = 1);
        // Write roots for given set of coefficients.
        for (size_t root_id = 0;
             root_id < conf_roots[coeff_root_id].size();
             ++root_id) {
          num_written = fwrite(
              reinterpret_cast<char const *>(
                  &conf_roots[coeff_root_id][root_id]),
              sizeof(conf_roots[coeff_root_id][root_id]), 1, fp);
          assert(num_written = 1);
        }
      }
    }
  }

  return 0;
}
