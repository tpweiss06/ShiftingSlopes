#!/bin/bash -l
#PBS -l walltime=01:00:00,nodes=2:ppn=24,pmem=2580mb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Change to the relevant working directory
cd ~/ShiftingSlopes/ShiftingRange

Rscript=4c_Extinctions.R
LogFile=4c_Extinctions.log

# Load R and MPI
module load R/3.4.4
module load ompi/3.0.0/gnu-7.2.0

export RMPI_TYPE=OPENMPI
export OMPI_MCA_mpi_warn_on_fork=0

# Run the R script itself, saving the output to
#	a log file
mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet $Rscript $LogFile
