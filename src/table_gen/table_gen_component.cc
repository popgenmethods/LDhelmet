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

#include "table_gen/table_gen_component.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <cstdlib>
#include <vector>

#include "common/read_confs.h"
#include "common/rho_finder.h"
#include "common/version_number.h"
#include "table_gen/output_writer.h"
#include "table_gen/solve_degree.h"
#include "table_gen/solve_degree-inl.h"
#include "table_gen/table_gen_options.h"

int TableGenComponent(std::string const base_command, int argc, char **argv) {
  uint64_t version_number = TABLE_GEN_VERSION_OUTPUT;
  uint64_t version_bit_string = TABLE_GEN_SALT + version_number;

  std::string version_string = TABLE_GEN_VERSION_STRING;

  CmdLineOptionsTableGen cmd_line_options(base_command,
                                          argc,
                                          argv,
                                          version_string);

  if (!cmd_line_options.success()) {
    std::exit(1);
  }

  printf("%s\n", ShowCommandLineOptions(argc, argv).c_str());

  uint32_t const num_threads = cmd_line_options.num_threads_;
  printf("Number of threads: %d.\n",
         static_cast<int>(num_threads));
  assert(num_threads > 0);

  std::string const &conf_file = cmd_line_options.conf_file_;
  std::string const &output_file = cmd_line_options.output_file_;

  std::vector<double> const &thetas = cmd_line_options.thetas_;

  assert(thetas.size() > 0);
  if (thetas.size() > 1) {
    fprintf(stderr,
            "Not implemented yet: Cannot handle more than one theta.\n");
    std::exit(1);
  }

  boost::tuple<std::vector<double>, std::vector<double> > rho_grid =
    ParseRhoRange(cmd_line_options.rho_range_);
  std::vector<double> rhos = GetRhoList(rho_grid);
  assert(rhos.size() > 0);

  if (rhos.front() != 0.0) {
    fprintf(stderr,
            "The rho values must begin at 0.0.\n");
    std::exit(1);
  }

  // Ensure that rhos is sorted.
  {
    std::vector<double> debug_rhos = rhos;
    sort(debug_rhos.begin(), debug_rhos.end());
    assert(rhos == debug_rhos);
  }

  assert(rhos.size() > 0);
  if (rhos.front() != 0.0) {
    fprintf(stderr,
            "The range of rho values must begin with 0.0.\n");
    std::exit(1);
  }

  printf("theta: ");
  for (int i = 0; i < static_cast<int>(thetas.size()); ++i) {
    printf("%g ", thetas[i]);
  }
  printf("\n");

  printf("rho values: ");
  for (int i = 0; i < static_cast<int>(rhos.size()); ++i) {
    printf("%g ", rhos[i]);
  }
  printf("\n");

  // Load confs from file.

  printf("Reading configurations from file.\n");

  std::vector<Conf> const conf_list = LoadConfigurations(conf_file);

  std::vector<size_t> degree_seps;
  uint32_t max_degree, max_sample_size, max_locus;
  boost::tie(degree_seps, max_degree, max_sample_size, max_locus) =
    PreprocessConfs(conf_list);

  printf("Largest sample size: %d.\n", static_cast<int>(max_sample_size));

  {
    // Write input confs.
    InputConfBinaryWriter input_conf_binary_writer(output_file,
                                                   conf_list,
                                                   degree_seps);
    {
      // Write version number.
      int num_written;
      num_written = fwrite(reinterpret_cast<void const *>(&version_bit_string),
                           sizeof(version_bit_string),
                           1, input_conf_binary_writer.fp_);
      assert(num_written == 1);

      // Write num confs.
      uint64_t num_conf_list = conf_list.size();
      num_written = fwrite(reinterpret_cast<void const *>(&num_conf_list),
                           sizeof(num_conf_list),
                           1, input_conf_binary_writer.fp_);
      assert(num_written == 1);

      // Write theta.
      double theta = thetas.front();
      num_written = fwrite(reinterpret_cast<void const *>(&theta),
                           sizeof(theta),
                           1, input_conf_binary_writer.fp_);
      assert(num_written == 1);

      // If using LDHelmet to generate lookup tables turn off interpolation
      uint8_t interpolate = 0;
      num_written = fwrite(reinterpret_cast<void const *>(&interpolate),
                           sizeof(interpolate),
                           1, input_conf_binary_writer.fp_);
      assert(num_written == 1);

      // Write number of rho segments, excluding end point.
      if (rho_grid.get<0>().size() == 0) {
        fprintf(stderr,
                "Error: The size of rho_grid is 0. "
                "The error should have been caught earlier. "
                "There is probably a bug in the code.\n");
        std::exit(1);
      }
      uint64_t num_rho_segments = rho_grid.get<0>().size() - 1;
      num_written = fwrite(reinterpret_cast<void const *>(&num_rho_segments),
                           sizeof(num_rho_segments),
                           1, input_conf_binary_writer.fp_);
      assert(num_written == 1);

      // Write rho segments.
      assert(rho_grid.get<0>().size() > 0);
      for (uint32_t i = 0; i < rho_grid.get<0>().size() - 1; ++i) {
        double start = rho_grid.get<0>()[i];
        double delta = rho_grid.get<1>()[i];
        num_written = fwrite(reinterpret_cast<void const *>(&start),
                             sizeof(start),
                             1, input_conf_binary_writer.fp_);
        assert(num_written == 1);

        num_written = fwrite(reinterpret_cast<void const *>(&delta),
                             sizeof(delta),
                             1, input_conf_binary_writer.fp_);
        assert(num_written == 1);
      }
      num_written = fwrite(reinterpret_cast<void const *>(
                               &(rho_grid.get<0>().back())),
                           sizeof(rho_grid.get<0>().back()),
                           1, input_conf_binary_writer.fp_);
      assert(num_written == 1);

      // Write confs.
      for (std::vector<Conf>::const_iterator citer = conf_list.begin();
           citer != conf_list.end();
           ++citer) {
        Conf::BinaryRep::RepType rep = Conf::BinaryRep(*citer).rep();
        num_written = fwrite(reinterpret_cast<void const *>(&rep),
                             sizeof(rep),
                             1, input_conf_binary_writer.fp_);
        assert(num_written == 1);
      }
    }

    SCCMutation const SCCPred(max_locus);

    assert(thetas.size() == 1);
    for (size_t theta_id = 0; theta_id < thetas.size(); ++theta_id) {
      for (size_t rho_id = 0; rho_id < rhos.size(); ++rho_id) {
        printf("rho: %g\n", rhos[rho_id]);

        Vec8 table(max_degree + 1);

        printf("Degrees: ");

        for (uint64_t degree = 2; degree <= max_degree; ++degree) {
            // Solve given degree.
            SolveDegree(thetas[theta_id],
                        rhos[rho_id],
                        &table,
                        degree,
                        max_degree,
                        SCCPred,
                        num_threads);

            printf("%d ", static_cast<int>(degree));
            fflush(stdout);

            // Write out likelihoods.
            input_conf_binary_writer.Write(degree, table);
        }
        printf("\n");
      }
    }
  }

  return 0;
}
