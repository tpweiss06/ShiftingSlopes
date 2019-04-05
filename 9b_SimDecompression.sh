#!/bin/bash -l
#PBS -N FolderDecompression
#PBS -l walltime=06:00:00,nodes=1:ppn=1,mem=2580mb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Change to the relevant working directory
cd ~/ShiftingSlopes


# zip all the simulation folders once the data extraction scripts have run
gunzip -r MainSim/
gunzip -r Slow/
gunzip -r Fast/

