# This script will extract the extinction information from the simulations

setwd("~/ShiftingSlopes/ShiftingRange")
library(parallel)
library(Rmpi)
nProc <- 24*2
source("ShiftParams.R")

# Calculate the time period to assess extinctions over
TotalTime <- BurnIn + LengthShift + BurnOut

# Create an object to hold the extinction data
Extinctions <- array(NA, dim = c(3, 9, TotalTime))

SimVec <- NULL
for(v in 1:3){
     for(p in 1:9){
          FilePath <- paste("~/ShiftingSlopes/ShiftingRange/", SpeedWords[v],
                            "/Params", p, sep = "")
          SimFiles <- list.files(FilePath)
          FullPaths <- paste(FilePath, "/", SimFiles, "/SummaryStats.csv", sep = "")
          SimVec <- c(SimVec, FullPaths)
     }
}

# Create a function to extract the abundance information I want from each simulation
ExtinctExtract <- function(InFile){
     SimData <- read.csv(InFile)
     GenExtinct <- max(SimData$gen) + 1    # Generations stop being recorded at 0 population
     return(GenExtinct)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")
ExtGens <- clusterApply(cl, x = SimVec, fun = ExtinctExtract)

for(v in 1:3){
     for(p in 1:9){
          CurPattern <- paste(SpeedWords[v], "/Params", p, sep = "")
          CurIndices <- grep(pattern = CurPattern, x = SimVec)
          CurExtGens <- rep(NA, length(CurIndices))
          for(i in 1:length(CurIndices)){
               CurExtGens[i] <- ExtGens[[CurIndices[i]]]
          }
          for(t in 1:TotalTime){
               Extinctions[v,p,t] <- sum(CurExtGens == t)
          }
     }
}

# Save the extinction output
save(Extinctions, file = "Extinctions.rdata")

# The code will then shut down the cluster after this line.
