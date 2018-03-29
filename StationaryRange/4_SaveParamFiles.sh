#!/bin/bash -l
#PBS -l walltime=00:00:05,nodes=1:ppn=1,pmem=2580mb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Params1
#cd ~/ShiftingSlopes/StationaryRange/Params1
#SubFolders=($( ls ))
#cp ${SubFolders[1]}/parameters.R > ~/ShiftingSlopes/ShiftingRange/ParamFiles/Params1.R

# Params2
#cd ~/ShiftingSlopes/StationaryRange/Params2
#SubFolders=($( ls ))
#cp ${SubFolders[1]}/parameters.R > ~/ShiftingSlopes/ShiftingRange/ParamFiles/Params2.R

# Params3
#cd ~/ShiftingSlopes/StationaryRange/Params3
#SubFolders=($( ls ))
#cp ${SubFolders[1]}/parameters.R > ~/ShiftingSlopes/ShiftingRange/ParamFiles/Params3.R

for (( i=1; i<10; i++));
do
cd ~/ShiftingSlopes/StationaryRange/Params$i
SubFolders=($( ls ))
cp ${SubFolders[1]}/parameters.R > ~/ShiftingSlopes/ShiftingRange/ParamFiles/Params${i}.R
done
