#!/bin/bash -l
#PBS -l walltime=00:45:00,nodes=1:ppn=1,mem=2gb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Set the names of R scripts and log files
Rscript=1a_SimCheck.R
LogFile=1a_SimCheck.log

# Change to the relevant working directory
cd ~/ShiftingSlopes/

# Load R and MPI
module load R/3.4.4
module load ompi/3.0.0/gnu-7.2.0-centos7

export RMPI_TYPE=OPENMPI
export OMPI_MCA_mpi_warn_on_fork=0

mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet $Rscript $LogFile
