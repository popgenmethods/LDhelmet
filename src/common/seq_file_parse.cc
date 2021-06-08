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

#include "common/seq_file_parse.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>

int const kBufferSize = 1024;

int GetLine(FILE *fp, std::string *string_buffer) {
  *string_buffer = std::string();
  bool not_eof = true;
  char char_buffer[kBufferSize];
  while (!feof(fp)) {
    char *fgets_ret = fgets(char_buffer, kBufferSize, fp);
    if (!feof(fp)) {
      assert(!ferror(fp));
      int char_buffer_len = strlen(char_buffer);
      assert(char_buffer_len > 0);
      if (char_buffer[char_buffer_len - 1] == '\n') {
        char_buffer[char_buffer_len - 1] = '\0';  // Remove newline.
        *string_buffer += char_buffer;
        break;
      } else {
        *string_buffer += char_buffer;
      }
    } else {
      assert(fgets_ret == NULL);
      not_eof = false;
    }
  }
  assert(!ferror(fp));
  return not_eof;
}

std::vector<std::string> ReadFastaFile(FILE *fp) {
  std::vector<std::string> seqs;

  std::string buffer;
  while (GetLine(fp, &buffer)) {
    if (buffer[0] == '>') {
      seqs.push_back(std::string());
    } else {
      seqs.back() += buffer;
    }
  }

  if (ferror(fp)) {
    fprintf(stderr,
            "Error reading input.\n");
    std::exit(1);
  }

  // Check if final line of FASTA file is a comment line,
  // which would be bad formatting.
  if (buffer[0] == '>') {
    fprintf(stderr,
            "File not in proper FASTA format.\n");
    std::exit(1);
  }

  return seqs;
}

std::vector<std::string> ReadFastqFile(FILE *fp) {
  std::vector<std::string> seqs;

  std::string buffer0, buffer1, buffer2, buffer3;

  while (true) {
    GetLine(fp, &buffer0);
    bool line_success0 = !feof(fp) && !ferror(fp);
    GetLine(fp, &buffer1);
    bool line_success1 = !feof(fp) && !ferror(fp);
    GetLine(fp, &buffer2);
    bool line_success2 = !feof(fp) && !ferror(fp);
    GetLine(fp, &buffer3);
    bool line_success3 = !feof(fp) && !ferror(fp);

    // Termination condition is for all the lines to fail reading.
    // There also should be no reading error.
    if (!ferror(fp)    &&
        !line_success0 &&
        !line_success1 &&
        !line_success2 &&
        !line_success3) {
      break;
    }

    if (!(line_success0 &&
          line_success1 &&
          line_success2 &&
          line_success3)
        || buffer0[0] != '@'
        || buffer2[0] != '+') {
      fprintf(stderr, "Input not in proper FASTQ format.\n");
      std::exit(1);
    }

    if (ferror(fp)) {
      fprintf(stderr, "Error in reading FASTQ file.\n");
      std::exit(1);
    }

    seqs.push_back(buffer1);
  }

  if (ferror(fp)) {
    fprintf(stderr, "Error reading input.\n");
    std::exit(1);
  }

  return seqs;
}

std::vector<std::string> ReadSeqFile(std::string const &input_file) {
  FILE *fp = fopen(input_file.c_str(), "r");
  if (fp == NULL) {
    fprintf(stderr, "Could not open file %s.\n", input_file.c_str());
    std::exit(1);
  }

  enum {FASTA, FASTQ} file_type;

  // Check first line to determine if its FASTA or FASTQ.
  {
    std::string buffer;
    GetLine(fp, &buffer);
    if (feof(fp) || ferror(fp)) {
      fprintf(stderr, "Error reading input when checking file type.\n");
      std::exit(1);
    }
    if (buffer[0] == '>') {
      file_type = FASTA;
    } else if (buffer[0] == '@') {
      file_type = FASTQ;
    } else {
      fprintf(stderr,
              "Input not in proper FASTA or FASTQ format.\n");
      std::exit(1);
    }

    rewind(fp);
  }

  if (file_type == FASTA) {
    return ReadFastaFile(fp);
  } else if (file_type == FASTQ) {
    return ReadFastqFile(fp);
  } else {
    fprintf(stderr,
            "There is a bug in the sequence file parser.\n");
    std::exit(1);
  }

  assert(fp != NULL);
  fclose(fp);
}

std::vector<uint64_t> ReadPosFile(std::string const &pos_file) {
  FILE *fp = fopen(pos_file.c_str(), "r");
  if (fp == NULL) {
    fprintf(stderr,
            "Could not open SNP positions file: %s.\n", pos_file.c_str());
    std::exit(1);
  }

  std::vector<uint64_t> snp_pos;
  uint64_t snp_pos1;

  int tmp_snp_pos1;
  int number_read = fscanf(fp, "%d", &tmp_snp_pos1);
  assert(number_read == 1);
  snp_pos1 = tmp_snp_pos1;
  snp_pos.push_back(snp_pos1);


  while (!feof(fp)) {
    number_read = fscanf(fp, "%d", &tmp_snp_pos1);
    if (number_read == 1) {
      snp_pos1 = tmp_snp_pos1;
      if (snp_pos1 <= snp_pos.back()) {
        fprintf(stderr,
                "Error: SNP positions must be in strictly increasing order "
                "and non-negative.\n");
        std::exit(1);
      }
      snp_pos.push_back(snp_pos1);
    } else {
      assert(number_read == -1 && feof(fp));
    }
  }

  assert(fp != NULL);
  fclose(fp);
  return snp_pos;
}
