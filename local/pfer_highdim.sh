#!/bin/bash
AMP=$1
SEED=$2

ml R/3.5

Rscript --vanilla ../simulations/pfer_highdim.R $AMP $SEED

