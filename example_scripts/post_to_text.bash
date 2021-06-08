#!/bin/bash

mkdir -p output &&
time ../ldhelmet post_to_text -m -p 0.025 -p 0.50 -p 0.975 -o output/output.txt output/output.post
