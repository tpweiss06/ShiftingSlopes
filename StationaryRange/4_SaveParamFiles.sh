#!/bin/bash -l
#PBS -l walltime=00:00:05,nodes=1:ppn=1,pmem=2580mb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

for (( i=1; i<10; i++));
do
cd ~/ShiftingSlopes/StationaryRange/Params$i
SubFolders=($( ls ))
cp ${SubFolders[1]}/parameters.R ~/ShiftingSlopes/ShiftingRange/ParamFiles/Params${i}.R
done
