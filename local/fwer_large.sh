#!/bin/bash
AMP=$1
SEED=$2

ml R/3.5

Rscript --vanilla ../simulations/fwer_large.R $AMP $SEED

