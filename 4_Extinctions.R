# This script will extract the extinction information from the simulations

setwd("~/ShiftingSlopes/")
library(parallel)
library(Rmpi)
nProc <- 24*2

# Calculate the time period to assess extinctions over
LengthShift <- 100
TotalTime <- 2150
BurnIn <- 2000
NumSims <- 200

# Create two objects to hold the extinction data: 1) an array to hold the number
#    of extinctions through time to use for the plots of the probability of extinction
#    and 2) a data frame with each simulation ID and a column of 1's and 0's for
#    each speed of climate change with 1 indicating extant and 0 indicating extinct
Extinctions <- array(NA, dim = c(3, 9, LengthShift))
SimIDs <- data.frame(SimID = rep(NA, NumSims*9), Params = rep(NA, NumSims*9),
                     Slow = rep(NA, NumSims*9), Moderate = rep(NA, NumSims*9), 
                     Fast = rep(NA, NumSims*9))

# First, populate the SimID column of the above data frame
k <- 1
for(p in 1:9){
     FilePath <- paste("~/ShiftingSlopes/MainSim/Params", p, sep = "")
     SimFiles <- list.files(FilePath)
     for(i in 1:NumSims){
          CurSimID <- strsplit(x = SimFiles[i], split = "/")[[1]][5]
          SimIDs$SimID[k] <- CurSimID
          SimIDs$Params[k] <- p
          k <- K + 1
     }
}

SimVec <- NULL
SpeedWords <- c("Slow", "MainSim", "Fast")
for(v in 1:3){
     for(p in 1:9){
          FilePath <- paste("~/ShiftingSlopes/", SpeedWords[v], "/Params", p, sep = "")
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
          if(v == 2){
               CurExtGens <- CurExtGens - BurnIn
          }
          for(t in 1:LengthShift){
               Extinctions[v,p,t] <- sum(CurExtGens == t)
          }
     }
}

# Save the extinction output
save(Extinctions, file = "Extinctions.rdata")

# Also, use the function output to construct a data frame of simulation IDs,
#    parameter combinations, and columns for slow, moderate, and fast speeds of
#    climate change with 1 for extant and 0 for extinct
for(i in 1:nrow(SimIDs)){
     SlowPattern <- paste("Slow/Params", SimIDs$Params[i], "/", SimIDs$SimID[i], sep = "")
     ModeratePattern <- paste("MainSim/Params", SimIDs$Params[i], "/", SimIDs$SimID[i], sep = "")
     FastPattern <- paste("Fast/Params", SimIDs$Params[i], "/", SimIDs$SimID[i], sep = "")
     SlowIndex <- grep(pattern = SlowPattern, x = SimVec)
     ModerateIndex <- grep(pattern = ModeratePattern, x = SimVec)
     FastIndex <- grep(pattern = FastPattern, x = SimVec)
     if(ExtGens[[SlowIndex]] >= LengthShift){
          SimIDs$Slow[i] <- 1
     }else{
          SimIDs$Slow[i] <- 0
     }
     if(ExtGens[[ModerateIndex]] >= (BurnIn + LengthShift) ){
          SimIDs$Moderate[i] <- 1
     }else{
          SimIDs$Moderate[i] <- 0
     }
     if(ExtGens[[FastIndex]] >= LengthShift){
          SimIDs$Fast[i] <- 1
     }else{
          SimIDs$Fast[i] <- 0
     }
}

save(SimIDs, file = "ExtSimIDs.rdata")

# The code will then shut down the cluster after this line.
