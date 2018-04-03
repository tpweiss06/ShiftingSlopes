#!/bin/bash -l
#PBS -l walltime=09:00:00,nodes=1:ppn=1,pmem=2580mb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

for (( i=1; i<10; i++));
do
     cd ~/ShiftingSlopes/StationaryRange/Params$i
     SubFolders=($( ls ))
     for (( j=1; j<101; j++));
     do
          gzip ${SubFolders[$j]}/PopMat.csv
     done
done

