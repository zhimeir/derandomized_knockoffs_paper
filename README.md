# derandomized_knockoffs_paper
A repository to reproduce the numerical results in the [Derandomized Knockoffs]() paper. 

## Overview
We provide the code to <em>exactly</em> reproduce the numerical results in Section 3-5 of the paper. Based on the  code,
we additionally develop an R package `derandomKnock`, which will be constantly improved and udpated. We recommend the users to check
out the R package unless their aim is to reproduce the examples of our paper.

## Folders
- `R/`: contains the main functions that implement the variable selection procedures.
- `simulations/`: contains the scripts to carry out the simulations.
- `results/`: stores the result.
- `bash/`: bash files to run the simulations.

## Usage
### One run
Each script in the `simulations/` folder implements one run of the simulation. The users can specify the amplitude and the seed when running the script. 
For example, to implement one run of the small-scale experiment in Section 3 with amplitude 4 and seed 1, run the following command in your terminal:
```{r}
Rscript ./simulations/pfer_small.R 4 1
```

### Multiple runs
The results presented in the paper are averaged over multiple runs. The user can use the bash file in `bash/` to automatically
implement multiple runs of the simulation. But note that it may take a long time if it is run on a laptop. To use the bash file,
run the following code in  your terminal:
```{r}
bash ./bash/run_all.sh
```
