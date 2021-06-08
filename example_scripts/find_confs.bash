#!/bin/bash

mkdir -p output &&
../ldhelmet find_confs --num_threads 24 -w 50 -o output/output.conf input.fasta
