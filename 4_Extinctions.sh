#!/bin/bash -l
#PBS -l walltime=00:30:00,nodes=2:ppn=24,mem=50gb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Change to the relevant working directory
cd ~/ShiftingSlopes/

# Set the names of R scripts and log files
Rscript=4_Extinctions.R
LogFile=4_Extinctions.log

# Load R and MPI
module load R/3.4.4
module load ompi/3.0.0/gnu-7.2.0-centos7

export RMPI_TYPE=OPENMPI
export OMPI_MCA_mpi_warn_on_fork=0

# Run the R script itself, saving the output to
#	a log file
mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet $Rscript $LogFile
