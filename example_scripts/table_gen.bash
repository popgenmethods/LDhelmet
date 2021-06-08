#!/bin/bash

mkdir -p output &&
time ../ldhelmet table_gen --num_threads 24 -t 0.008 -r 0.0 0.1 10.0 1.0 100.0 -c output/output.conf -o output/output.lk 
