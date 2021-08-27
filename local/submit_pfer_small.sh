#!/bin/bash
AMP_LIST=("3" "4" "5" "6" "7" "8")

# Slurm parameters
PART=candes,hns,stat,pilanci      # Partition names
MEMO=1G                         # Memory required (1GB)
TIME=00-02:00:00                  # Time required (2HRS)
CORE=10                           # Cores required (10)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
RES_DIR="../results/pfer_small"
mkdir -p $RES_DIR
LOGS="logs"
mkdir -p $LOGS

# Loop through the configurations
for AMP in "${AMP_LIST[@]}"; do
  for SEED in {1..100}; do

    # Script to be run
    SCRIPT="pfer_small.sh $AMP $SEED"

    # Define job name for this chromosome
    JOBN="pfer_small_"$AMP"_"$SEED
    OUTF=$LOGS"/"$JOBN".out"
    ERRF=$LOGS"/"$JOBN".err"

    # Assemble slurm order for this job
    ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT

    #Print order
    echo $ORD

    # Submit order
    $ORD

  done
done




