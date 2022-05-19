#!/bin/bash

if [ -z "$LD_LIBRARY_PATH" ]; then
  LD_LIBRARY_PATH="/gpfs0/biores/users/mishmarlab/Hadar/mymodule/gcc-7.5.0/lib64"
else
  LD_LIBRARY_PATH="/gpfs0/biores/users/mishmarlab/Hadar/mymodule/gcc-7.5.0/lib64:$LD_LIBRARY_PATH:/gpfs0/biores/apps/Miniconda2/Miniconda_v4.3.21/lib"
fi
export LD_LIBRARY_PATH 

/gpfs0/biores/users/mishmarlab/Anne/R-3.6.3/bin/Rscript  /gpfs0/biores/users/mishmarlab/Anne/1_DATASET/codi/STAT/bootstrap.r
