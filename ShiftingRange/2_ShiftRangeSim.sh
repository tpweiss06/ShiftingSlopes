#!/bin/bash -l
#PBS -l walltime=06:00:00,nodes=18:ppn=24,pmem=2580mb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Set the parameter combination to use and generate names of R scripts and log files
Rscript=2_ShiftRangeSim.R
LogFile=2_ShiftRangeSim.log

# Change to the relevant working directory
cd ~/ShiftingSlopes/ShiftingRange

# Load R and MPI
module load R/3.4.4
module load ompi/3.0.0/gnu-7.2.0

export RMPI_TYPE=OPENMPI
export OMPI_MCA_mpi_warn_on_fork=0

mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet $Rscript $LogFile
