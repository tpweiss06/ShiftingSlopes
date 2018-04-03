#!/bin/bash -l
#PBS -l walltime=00:05:00,nodes=1:ppn=1,pmem=2580mb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Change to the relevant working directory
cd ~/ShiftingSlopes/ShiftingRange

# Load R
module load R/3.4.4

# Run the R script itself, saving the output to
#	a log file
R --vanilla < SimCleanUp.R >& SimCleanUp.log
