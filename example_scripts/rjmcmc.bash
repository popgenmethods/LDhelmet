#!/bin/bash

mkdir -p output &&
time ../ldhelmet rjmcmc --num_threads 24 -l output/output.lk -p output/output.pade -m mut_mat.txt -s input.fasta -b 50 --burn_in 10000 -n 100000 -o output/output.post
