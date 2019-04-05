#!/bin/bash -l
#PBS -l walltime=03:00:00,nodes=4:ppn=24,mem=248gb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Set the names of R scripts and log files
Rscript=6_DispData.R
LogFile=6_DispData.log

# Change to the relevant working directory
cd ~/ShiftingSlopes/

# Load R and MPI
module load R/3.4.4
module load ompi/3.0.0/gnu-7.2.0-centos7

export RMPI_TYPE=OPENMPI
export OMPI_MCA_mpi_warn_on_fork=0

mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet $Rscript $LogFile
