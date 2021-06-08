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

#include "post_to_text/post_to_text_component.h"

#include <stdint.h>
#include <stdio.h>

#include <cassert>
#include <cstdlib>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include "common/version_number.h"
#include "post_to_text/post_to_text.h"
#include "post_to_text/post_to_text_options.h"

int PostToTextComponent(std::string const base_command,
                        int argc,
                        char **argv) {
  std::string version_string = POST_TO_TEXT_VERSION_STRING;
  CmdLineOptionsPostToText cmd_line_options(base_command,
                                            argc,
                                            argv,
                                            version_string);

  if (!cmd_line_options.success()) {
    std::exit(1);
  }

  printf("%s\n", ShowCommandLineOptions(argc, argv).c_str());

  bool const &mean_p = cmd_line_options.mean_p_;
  std::vector<double> const &percs = cmd_line_options.percs_;
  std::string const &output_file = cmd_line_options.output_file_;
  std::vector<std::string> const &input_files = cmd_line_options.input_files_;

  assert(input_files.size() == 1);
  FILE *fp_in = fopen(input_files.front().c_str(), "r");
  if (fp_in == NULL) {
    fprintf(stderr, "Error opening input file.\n");
    std::exit(1);
  }

  for (size_t perc_id = 0; perc_id < percs.size(); ++perc_id) {
    if (percs[perc_id] < 0.0 || percs[perc_id] > 1.0) {
      fprintf(stderr,
              "Percentiles must be between 0.0 and 1.0, inclusive.\n");
      std::exit(1);
    }
  }

  RjmcmcParams params;
  RjmcmcResults results;
  boost::tie(params, results) = ReadPostData(fp_in);

  int const perc_precision = 3;
  int const rate_precision = 4;

  FILE *fp_out = fopen(output_file.c_str(), "w");
  if (fp_out == NULL) {
    fprintf(stderr, "Error opening output file.\n");
    std::exit(1);
  }

  std::string header0 = "# parameters";
  std::string params_string = Show(params);

  fprintf(fp_out, "%s\n", header0.c_str());
  fprintf(fp_out, "%s\n", params_string.c_str());

  std::string header1 = "# left_snp right_snp";
  fprintf(fp_out, "%s", header1.c_str());

  if (mean_p) {
    fprintf(fp_out, "%s", " mean");
  }

  for (size_t perc_id = 0; perc_id < percs.size(); ++perc_id) {
    fprintf(fp_out, " p%.*f", perc_precision, percs[perc_id]);
  }

  fprintf(fp_out, "\n");

  // Process each partition separately.
  for (size_t part_id = 0; part_id < params.num_snp_partitions_; ++part_id) {
    // Read partition.
    uint64_t hash_value;  // Checksum.
    int num_read;
    num_read = fread(reinterpret_cast<char *>(&hash_value),
                    sizeof(hash_value), 1, fp_in);
    assert(num_read == 1);

    uint64_t num_samps_part;
    num_read = fread(reinterpret_cast<char *>(&num_samps_part),
                     sizeof(num_samps_part), 1, fp_in);
    assert(num_read == 1);
    assert(num_samps_part == params.num_samps_);

    uint64_t begin_part = params.snp_partitions_[part_id].first;
    uint64_t end_part   = params.snp_partitions_[part_id].second;

    uint64_t len_part = end_part - begin_part;
    assert(len_part >= 2);

    // [snpID, samp_id]
    std::vector<std::vector<double> > samps(
        len_part - 1,
        std::vector<double>(num_samps_part));

    for (size_t samp_id = 0; samp_id < num_samps_part; ++samp_id) {
      uint64_t num_cps;
      num_read = fread(reinterpret_cast<char *>(&num_cps),
                       sizeof(num_cps), 1, fp_in);
      assert(num_read == 1);
      assert(num_cps >= 2);

      std::vector<uint64_t> cp_pos(num_cps);
      std::vector<double> cp_rates(num_cps);
      for (size_t cp_id = 0; cp_id < num_cps; ++cp_id) {
        uint64_t snp_pos1;
        double rate;

        num_read = fread(reinterpret_cast<char *>(&snp_pos1),
                         sizeof(snp_pos1), 1, fp_in);
        assert(num_read == 1);
        num_read = fread(reinterpret_cast<char *>(&rate),
                         sizeof(rate), 1, fp_in);
        assert(num_read == 1);

        cp_pos[cp_id] = snp_pos1;
        cp_rates[cp_id] = rate;
      }

      assert(cp_pos.front() == params.snp_pos_[begin_part]);
      assert(cp_pos.back() == params.snp_pos_[end_part - 1]);

      uint64_t snp_idx = 0;
      for (size_t cp_id = 0; cp_id < num_cps - 1; ++cp_id) {
        assert(cp_id + 1 < cp_pos.size());
        assert(begin_part + snp_idx < params.snp_pos_.size());
        while (params.snp_pos_[begin_part + snp_idx] < cp_pos[cp_id + 1]) {
          assert(snp_idx < samps.size());
          assert(cp_id < cp_rates.size());
          samps[snp_idx][samp_id] = cp_rates[cp_id];
          ++snp_idx;
          assert(begin_part + snp_idx < params.snp_pos_.size());
        }
      }
    }

    // Write to file.
    for (size_t i = 0; i < len_part - 1; ++i) {
      fprintf(fp_out, "%d %d",
              static_cast<int>(params.snp_pos_[begin_part + i]),
              static_cast<int>(params.snp_pos_[begin_part + i + 1]));

      if (mean_p) {
        fprintf(fp_out, " %.*e",
                rate_precision, results.mean_[begin_part + i]);
      }

      // Write out percentiles.
      for (size_t perc_id = 0; perc_id < percs.size(); ++perc_id) {
        assert(i < samps.size());
        double perc = FindPerc(samps[i], percs[perc_id]);

        fprintf(fp_out, " %.*e",
                rate_precision, perc);
      }
      fprintf(fp_out, "\n");
    }
  }

  assert(fp_in != NULL);
  fclose(fp_in);

  assert(fp_out != NULL);
  fclose(fp_out);

  return 0;
}
