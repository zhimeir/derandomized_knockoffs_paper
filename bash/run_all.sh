#!/bin/bash
for amp in {3..8}
do
  for seed in {1..200}
  do
    Rscript ../simulations/pfer_small.R $amp $seed
    Rscript ../simulations/pfer_large.R $amp $seed
  done
done

amp_list=("6" "8" "10" "12" "14" "16")
for amp in ${amp_list[@]}
do
  for seed in {1..200}
  do
    Rscript ../simulations/pfer_small.R $amp $seed
    Rscript ../simulations/pfer_large.R $amp $seed
  done
done
