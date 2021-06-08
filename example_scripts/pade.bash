#!/bin/bash

mkdir -p output &&
../ldhelmet pade --num_threads 24 -t 0.008 -x 12 --defect_threshold 40 -c output/output.conf -o output/output.pade
