#!/bin/bash

mkdir -p output &&
../ldhelmet max_lk --num_threads 24 -l output/output.lk -p output/output.pade -m mut_mat.txt -s input.fasta
