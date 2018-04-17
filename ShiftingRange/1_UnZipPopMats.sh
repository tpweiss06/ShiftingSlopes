#!/bin/bash -l
#PBS -l walltime=03:00:00,nodes=1:ppn=1,mem=2gb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

for (( i=1; i<10; i++));
do
     cd ~/ShiftingSlopes/StationaryRange/Params$i
     SubFolders=($( ls ))
     for (( j=1; j<101; j++));
     do
          gunzip ~/ShiftingSlopes/StationaryRange/Params${i}/${SubFolders[$j]}/PopMat.csv.gz
     done
done

