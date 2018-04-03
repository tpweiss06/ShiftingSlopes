#!/bin/bash -l
#PBS -N FolderCompression
#PBS -l walltime=18:00:00,nodes=1:ppn=1
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Change to the relevant working directory
cd ~/ShiftingSlopes/ShiftingRange/Slow

# Now sequentially gzip all the files in the simulation folders to save space
gzip -r Params1/
gzip -r Params2/
gzip -r Params3/
gzip -r Params4/
gzip -r Params5/
gzip -r Params6/
gzip -r Params7/
gzip -r Params8/
gzip -r Params9/

# Change to the relevant working directory
cd ~/ShiftingSlopes/ShiftingRange/Med

# Now sequentially gzip all the files in the simulation folders to save space
gzip -r Params1/
gzip -r Params2/
gzip -r Params3/
gzip -r Params4/
gzip -r Params5/
gzip -r Params6/
gzip -r Params7/
gzip -r Params8/
gzip -r Params9/

# Change to the relevant working directory
cd ~/ShiftingSlopes/ShiftingRange/Fast

# Now sequentially gzip all the files in the simulation folders to save space
gzip -r Params1/
gzip -r Params2/
gzip -r Params3/
gzip -r Params4/
gzip -r Params5/
gzip -r Params6/
gzip -r Params7/
gzip -r Params8/
gzip -r Params9/
