#!/bin/bash -l
#PBS -N FolderCompression
#PBS -l walltime=06:00:00,nodes=1:ppn=1,mem=2580mb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Change to the relevant working directory
cd ~/ShiftingSlopes


# zip all the simulation folders once the data extraction scripts have run
gzip -r MainSim/
gzip -r Slow/
gzip -r Fast/

