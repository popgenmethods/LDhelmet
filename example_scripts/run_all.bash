#!/bin/bash

mkdir -p output &&
bash find_confs.bash   &&
bash table_gen.bash    &&
bash pade.bash         &&
bash rjmcmc.bash       &&
bash post_to_text.bash &&
bash max_lk.bash
