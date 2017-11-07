#!/bin/bash -l
#PBS -N ExtractAbunds
#PBS -l walltime=00:20:00,nodes=1:ppn=24
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Change to the relevant working directory
cd ~/ShiftingSlopes/StationaryRange

# Load R and MPI
module load R/3.3.2
module load ompi/2.0.1/intel-2016-update3

export RMPI_TYPE=OPENMPI

# Run the R script itself, saving the output to
#	a log file
R --vanilla < ExtractAbunds.R >& ExtractAbunds.log 