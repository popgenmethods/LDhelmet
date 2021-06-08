# Copyright (C) 2012  Andrew H. Chan
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# <http://www.gnu.org/licenses/>
# 

########
# Change INC_FLAG and LIB_FLAG to point to your installations of the Boost C++
# Libraries and the GNU Scientific Library, if you performed a manual
# installation of these packages.
#
# Otherwise blank INC_FLAG and LIB_FLAG variables should suffice.
#
# Uncomment and change these variables as necessary.
#

INC_FLAG =
LIB_FLAG =
#INC_FLAG = -I/home/gradstud/myBoost/include -I/home/gradstud/myGSL/include
#LIB_FLAG = -L/home/gradstud/myBoost/lib -L/home/gradstud/myGSL/lib
LIB = -lboost_thread -lboost_program_options -lboost_system -lgsl -lgslcblas


# Typical for OS X:
#INC_FLAG = -I/usr/local/include/boost/ -I/usr/local/include/gsl/ 
#LIB_FLAG = -L/usr/local/lib
#LIB = -lboost_thread-mt -lboost_program_options-mt -lboost_system-mt -lgsl -lgslcblas


#
########

CC = g++
CFLAGS = -std=c++11 -O3 -Isrc $(INC_FLAG)
LFLAGS = $(LIB_FLAG)

all: build_dir ldhelmet

build_dir:
	mkdir -p build

# common

common: common_dir build/common/command_line_options.o build/common/conf.o build/common/mar.o build/common/snp_partitions.o build/common/mut_mat_prior.o build/common/seq_file_parse.o build/common/seq_process.o build/common/load_data.o build/common/conf_gen.o build/common/rho_finder.o build/common/lk_pade_table.o build/common/site_map_log_lk.o build/common/log_lk_computer.o build/common/ncr.o build/common/read_confs.o

common_dir:
	mkdir -p build/common

build/common/command_line_options.o: src/common/command_line_options.cc src/common/command_line_options.h
	$(CC) $(CFLAGS) -c -o build/common/command_line_options.o src/common/command_line_options.cc

build/common/conf.o: src/common/conf.cc src/common/conf.h
	$(CC) $(CFLAGS) -c -o build/common/conf.o src/common/conf.cc

build/common/mar.o: src/common/mar.cc src/common/mar.h src/common/conf.h
	$(CC) $(CFLAGS) -c -o build/common/mar.o src/common/mar.cc

build/common/snp_partitions.o: src/common/snp_partitions.cc src/common/snp_partitions.h
	$(CC) $(CFLAGS) -c -o build/common/snp_partitions.o src/common/snp_partitions.cc

build/common/mut_mat_prior.o: src/common/mut_mat_prior.cc src/common/mut_mat_prior.h
	$(CC) $(CFLAGS) -c -o build/common/mut_mat_prior.o src/common/mut_mat_prior.cc

build/common/seq_file_parse.o: src/common/seq_file_parse.cc src/common/seq_file_parse.h
	$(CC) $(CFLAGS) -c -o build/common/seq_file_parse.o src/common/seq_file_parse.cc

build/common/seq_process.o: src/common/seq_process.cc src/common/seq_process.h
	$(CC) $(CFLAGS) -c -o build/common/seq_process.o src/common/seq_process.cc

build/common/load_data.o: src/common/load_data.cc src/common/load_data.h src/common/mut_mat_prior.h src/common/seq_file_parse.h src/common/seq_process.h src/find_confs/find_confs.h
	$(CC) $(CFLAGS) -c -o build/common/load_data.o src/common/load_data.cc

build/common/conf_gen.o: src/common/conf_gen.cc src/common/conf_gen.h src/common/conf.h src/common/mar.h src/common/predicates.h
	$(CC) $(CFLAGS) -c -o build/common/conf_gen.o src/common/conf_gen.cc

build/common/rho_finder.o: src/common/rho_finder.cc src/common/rho_finder.h
	$(CC) $(CFLAGS) -c -o build/common/rho_finder.o src/common/rho_finder.cc

build/common/lk_pade_table.o: src/common/lk_pade_table.cc src/common/lk_pade_table.h src/common/version_number.h src/common/conf.h src/common/rho_finder.h src/find_confs/find_confs.h
	$(CC) $(CFLAGS) -c -o build/common/lk_pade_table.o src/common/lk_pade_table.cc

build/common/site_map_log_lk.o: src/common/site_map_log_lk.cc src/common/site_map_log_lk.h src/common/rho_finder.h src/common/lk_pade_table.h src/common/seq_process.h src/common/mut_mat_prior.h
	$(CC) $(CFLAGS) -c -o build/common/site_map_log_lk.o src/common/site_map_log_lk.cc

build/common/log_lk_computer.o: src/common/log_lk_computer.cc src/common/log_lk_computer.h src/common/change_point.h src/common/site_map_log_lk.h
	$(CC) $(CFLAGS) -c -o build/common/log_lk_computer.o src/common/log_lk_computer.cc

build/common/ncr.o: src/common/ncr.cc src/common/ncr.h
	$(CC) $(CFLAGS) -c -o build/common/ncr.o src/common/ncr.cc

build/common/read_confs.o: src/common/read_confs.cc src/common/read_confs.h src/common/conf.h
	$(CC) $(CFLAGS) -c -o build/common/read_confs.o src/common/read_confs.cc

# find_confs

find_confs: find_confs_dir build/find_confs/find_confs_component.o build/find_confs/find_confs_options.o build/find_confs/find_confs.o

find_confs_dir:
	mkdir -p build/find_confs

build/find_confs/find_confs_component.o: src/find_confs/find_confs_component.cc src/find_confs/find_confs_component.h
	$(CC) $(CFLAGS) -c -o build/find_confs/find_confs_component.o src/find_confs/find_confs_component.cc

build/find_confs/find_confs_options.o: src/find_confs/find_confs_options.cc src/find_confs/find_confs_options.h
	$(CC) $(CFLAGS) -c -o build/find_confs/find_confs_options.o src/find_confs/find_confs_options.cc

build/find_confs/find_confs.o: src/find_confs/find_confs.cc src/find_confs/find_confs.h
	$(CC) $(CFLAGS) -c -o build/find_confs/find_confs.o src/find_confs/find_confs.cc

# max_lk

max_lk: max_lk_dir build/max_lk/max_lk_component.o build/max_lk/max_lk_options.o

max_lk_dir:
	mkdir -p build/max_lk

build/max_lk/max_lk_component.o: src/max_lk/max_lk_component.cc src/max_lk/max_lk_component.h src/common/lk_pade_table.h src/common/load_data.h src/common/log_lk_computer.h src/common/mut_mat_prior.h src/common/site_map_log_lk.h src/common/version_number.h src/max_lk/max_lk_options.h
	$(CC) $(CFLAGS) -c -o build/max_lk/max_lk_component.o src/max_lk/max_lk_component.cc

build/max_lk/max_lk_options.o: src/max_lk/max_lk_options.cc src/max_lk/max_lk_options.h src/common/command_line_options.h
	$(CC) $(CFLAGS) -c -o build/max_lk/max_lk_options.o src/max_lk/max_lk_options.cc

# post_to_text

post_to_text: post_to_text_dir build/post_to_text/post_to_text_component.o build/post_to_text/post_to_text.o build/post_to_text/post_to_text_options.o

post_to_text_dir:
	mkdir -p build/post_to_text

build/post_to_text/post_to_text_component.o: src/post_to_text/post_to_text_component.cc src/post_to_text/post_to_text_component.h src/common/version_number.h src/post_to_text/post_to_text.h src/post_to_text/post_to_text_options.h
	$(CC) $(CFLAGS) -c -o build/post_to_text/post_to_text_component.o src/post_to_text/post_to_text_component.cc

build/post_to_text/post_to_text.o: src/post_to_text/post_to_text.cc src/post_to_text/post_to_text.h src/common/conf.h src/common/version_number.h
	$(CC) $(CFLAGS) -c -o build/post_to_text/post_to_text.o src/post_to_text/post_to_text.cc

build/post_to_text/post_to_text_options.o: src/post_to_text/post_to_text_options.cc src/post_to_text/post_to_text_options.h src/common/command_line_options.h
	$(CC) $(CFLAGS) -c -o build/post_to_text/post_to_text_options.o src/post_to_text/post_to_text_options.cc

# table_gen

table_gen: table_gen_dir build/table_gen/table_gen_options.o build/table_gen/linear_solve.o build/table_gen/solve_degree.o build/table_gen/solve_scc.o build/table_gen/table_gen_component.o build/table_gen/output_writer.o build/table_gen/table_management.o build/table_gen/single_locus.o

table_gen_dir:
	mkdir -p build/table_gen

build/table_gen/table_gen_options.o: src/table_gen/table_gen_options.cc src/table_gen/table_gen_options.h src/common/command_line_options.h
	$(CC) $(CFLAGS) -c -o build/table_gen/table_gen_options.o src/table_gen/table_gen_options.cc

build/table_gen/linear_solve.o: src/table_gen/linear_solve.cc src/table_gen/linear_solve.h
	$(CC) $(CFLAGS) -c -o build/table_gen/linear_solve.o src/table_gen/linear_solve.cc

build/table_gen/solve_scc.o: src/table_gen/solve_scc.cc src/table_gen/solve_scc.h src/table_gen/linear_solve.h src/common/vector_definitions.h src/common/conf.h src/common/ncr.h src/table_gen/table_management.h
	$(CC) $(CFLAGS) -c -o build/table_gen/solve_scc.o src/table_gen/solve_scc.cc

build/table_gen/solve_degree.o: src/table_gen/solve_degree.cc src/table_gen/solve_degree.h src/common/conf.h src/common/mar.h src/common/predicates.h src/common/vector_definitions.h src/table_gen/single_locus.h src/table_gen/table_management.h
	$(CC) $(CFLAGS) -c -o build/table_gen/solve_degree.o src/table_gen/solve_degree.cc

build/table_gen/table_gen_component.o: src/table_gen/table_gen_component.cc src/table_gen/table_gen_component.h src/table_gen/degree_task.h src/table_gen/solve_degree.h src/table_gen/solve_degree-inl.h src/table_gen/table_gen_options.h src/common/conf_gen.h src/common/vector_definitions.h src/table_gen/output_writer.h src/common/predicates.h src/common/threading.h src/common/read_confs.h src/common/rho_finder.h src/common/version_number.h
	$(CC) $(CFLAGS) -c -o build/table_gen/table_gen_component.o src/table_gen/table_gen_component.cc

build/table_gen/output_writer.o: src/table_gen/output_writer.cc src/table_gen/output_writer.h src/common/vector_definitions.h src/common/conf.h src/table_gen/table_management.h
	$(CC) $(CFLAGS) -c -o build/table_gen/output_writer.o src/table_gen/output_writer.cc

build/table_gen/table_management.o: src/table_gen/table_management.cc src/table_gen/table_management.h src/common/vector_definitions.h src/common/conf.h
	$(CC) $(CFLAGS) -c -o build/table_gen/table_management.o src/table_gen/table_management.cc

build/table_gen/single_locus.o: src/table_gen/single_locus.cc src/table_gen/single_locus.h src/common/conf.h src/common/ncr.h
	$(CC) $(CFLAGS) -c -o build/table_gen/single_locus.o src/table_gen/single_locus.cc

# convert_table

convert_table: convert_table_dir build/convert_table/convert_table_options.o build/convert_table/convert_table_component.o

convert_table_dir:
	mkdir -p build/convert_table

build/convert_table/convert_table_options.o: src/convert_table/convert_table_options.cc src/convert_table/convert_table_options.h src/common/command_line_options.h
	$(CC) $(CFLAGS) -c -o build/convert_table/convert_table_options.o src/convert_table/convert_table_options.cc

build/convert_table/convert_table_component.o: src/convert_table/convert_table_component.cc src/convert_table/convert_table_component.h src/convert_table/convert_table_options.h src/common/conf_gen.h src/table_gen/output_writer.h src/common/predicates.h src/common/read_confs.h src/common/rho_finder.h src/common/version_number.h src/table_gen/table_management.h
	$(CC) $(CFLAGS) -c -o build/convert_table/convert_table_component.o src/convert_table/convert_table_component.cc

# pade

pade: pade_dir build/pade/pade_options.o build/pade/pade_component.o build/pade/subtable.o build/pade/coeff.o build/pade/compute_g.o build/pade/one_locus.o

pade_dir:
	mkdir -p build/pade

build/pade/pade_options.o: src/pade/pade_options.cc src/pade/pade_options.h src/common/command_line_options.h
	$(CC) $(CFLAGS) -c -o build/pade/pade_options.o src/pade/pade_options.cc

build/pade/pade_component.o: src/pade/pade_component.cc src/pade/pade_component.h src/common/conf_gen.h src/common/read_confs.h src/common/threading.h src/common/version_number.h src/pade/coeff.h src/pade/compute_g.h src/pade/g_task.h src/pade/memory_manager.h src/pade/memory_manager-inl.h src/pade/pade_options.h src/pade/q_task.h src/pade/subtable.h
	$(CC) $(CFLAGS) -c -o build/pade/pade_component.o src/pade/pade_component.cc

build/pade/subtable.o: src/pade/subtable.cc src/pade/subtable.h src/pade/subtable-inl.h src/common/conf_gen.h src/common/vector_definitions.h
	$(CC) $(CFLAGS) -c -o build/pade/subtable.o src/pade/subtable.cc

build/pade/coeff.o: src/pade/coeff.cc src/pade/coeff.h src/pade/subtable.h src/common/ncr.h
	$(CC) $(CFLAGS) -c -o build/pade/coeff.o src/pade/coeff.cc

build/pade/compute_g.o: src/pade/compute_g.cc src/pade/compute_g.h src/pade/one_locus.h src/pade/subtable.h src/common/conf.h
	$(CC) $(CFLAGS) -c -o build/pade/compute_g.o src/pade/compute_g.cc

build/pade/one_locus.o: src/pade/one_locus.cc src/pade/one_locus.h src/pade/one_locus.h
	$(CC) $(CFLAGS) -c -o build/pade/one_locus.o src/pade/one_locus.cc

# rjmcmc

rjmcmc: rjmcmc_dir build/rjmcmc/rjmcmc_options.o build/rjmcmc/rjmcmc_component.o build/rjmcmc/handle_output.o build/rjmcmc/r_task.o build/rjmcmc/acceptance_log.o build/rjmcmc/post_rho_map.o build/rjmcmc/rjmcmc.o build/rjmcmc/proposals.o

rjmcmc_dir:
	mkdir -p build/rjmcmc

build/rjmcmc/rjmcmc_options.o: src/rjmcmc/rjmcmc_options.cc src/rjmcmc/rjmcmc_options.h src/common/command_line_options.h src/rjmcmc/rjmcmc_options.h src/rjmcmc/ran_num_gen.h
	$(CC) $(CFLAGS) -c -o build/rjmcmc/rjmcmc_options.o src/rjmcmc/rjmcmc_options.cc

build/rjmcmc/rjmcmc_component.o: src/rjmcmc/rjmcmc_component.cc src/rjmcmc/rjmcmc_component.h src/common/load_data.h src/common/seq_file_parse.h src/common/seq_process.h src/common/version_number.h src/rjmcmc/handle_output.h src/rjmcmc/r_task.h src/rjmcmc/ran_num_gen.h src/rjmcmc/rjmcmc.h src/rjmcmc/rjmcmc_options.h
	$(CC) $(CFLAGS) -c -o build/rjmcmc/rjmcmc_component.o src/rjmcmc/rjmcmc_component.cc

build/rjmcmc/handle_output.o: src/rjmcmc/handle_output.cc src/rjmcmc/handle_output.h src/common/snp_partitions.h src/common/mut_mat_prior.h src/rjmcmc/rjmcmc_options.h
	$(CC) $(CFLAGS) -c -o build/rjmcmc/handle_output.o src/rjmcmc/handle_output.cc

build/rjmcmc/r_task.o: src/rjmcmc/r_task.cc src/rjmcmc/r_task.h src/common/lk_pade_table.h src/common/mut_mat_prior.h src/rjmcmc/acceptance_log.h src/rjmcmc/handle_output.h src/rjmcmc/priors.h src/rjmcmc/proposals.h src/rjmcmc/ran_num_gen.h src/rjmcmc/rjmcmc.h 
	$(CC) $(CFLAGS) -c -o build/rjmcmc/r_task.o src/rjmcmc/r_task.cc

build/rjmcmc/acceptance_log.o: src/rjmcmc/acceptance_log.cc src/rjmcmc/acceptance_log.h
	$(CC) $(CFLAGS) -c -o build/rjmcmc/acceptance_log.o src/rjmcmc/acceptance_log.cc

build/rjmcmc/post_rho_map.o: src/rjmcmc/post_rho_map.cc src/rjmcmc/post_rho_map.h
	$(CC) $(CFLAGS) -c -o build/rjmcmc/post_rho_map.o src/rjmcmc/post_rho_map.cc

build/rjmcmc/rjmcmc.o: src/rjmcmc/rjmcmc.cc src/rjmcmc/rjmcmc.h src/common/binary_search.h src/common/load_data.h src/common/log_lk_computer.h src/common/seq_process.h src/rjmcmc/acceptance_log.h src/rjmcmc/post_rho_map.h src/rjmcmc/priors.h src/rjmcmc/proposals.h src/rjmcmc/ran_num_gen.h
	$(CC) $(CFLAGS) -c -o build/rjmcmc/rjmcmc.o src/rjmcmc/rjmcmc.cc

build/rjmcmc/proposals.o: src/rjmcmc/proposals.cc src/rjmcmc/proposals.h
	$(CC) $(CFLAGS) -c -o build/rjmcmc/proposals.o src/rjmcmc/proposals.cc

# ldhelmet

ldhelmet: build/ldhelmet.o common find_confs table_gen post_to_text max_lk pade rjmcmc convert_table
	$(CC) $(CFLAGS) $(LFLAGS) -o ldhelmet build/ldhelmet.o build/*/*.o $(LIB)

build/ldhelmet.o: src/ldhelmet.cc src/common/version_number.h src/find_confs/find_confs_component.h src/table_gen/table_gen_component.h src/post_to_text/post_to_text_component.h src/max_lk/max_lk_component.h src/pade/pade_component.h src/convert_table/convert_table_component.h
	$(CC) $(CFLAGS) -c -o build/ldhelmet.o src/ldhelmet.cc

clean:
	rm -v ldhelmet build/*/*.o build/*.o; rmdir -v build/* build
