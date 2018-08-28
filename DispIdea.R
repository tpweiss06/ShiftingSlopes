# This script will extract the initial dispersal values of individuals immediately
#    prior to climate change (DispInit) and the change in mean dispersal values at
#    the patch level (along with the relative x coordinate of the patch) to assess
#    the role of evolution (DispEvol).

# Set the number of processors available for this scrip
nProc <- 4*24

# Set the working directory
setwd("~/ShiftingSlopes/")
source("~/ShiftingSlopes/SimFunctions.R")
library(parallel)
library(Rmpi)
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")
NumSims <- 200

# Set some initial parameter values to assess evolutionary change in dispersal
Threshold <- 10
BetaInit <- 0
BurnIn <- 2000
BurnOut <- 50
LengthShift <- 100
dmax <- 1000
eta <- RangeParams$eta[1]
rho <- 0.5
ClimSpeeds <- c(0.5, 1, 2)
MaxGen <- c(151, 2151, 151)
BetaCoord <- matrix(NA, nrow = 3, ncol = BurnIn + LengthShift + BurnOut)
for(i in 1:3){
     BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
               eta = eta, v = ClimSpeeds[i]) / eta
     BetaCoord[i,] <- c(rep(BetaInit, BurnIn), BetaShift, rep(BetaShift[LengthShift], BurnOut))
}

# Create a data frame to index the individual simulations and parameter values
TraitIndices <- expand.grid(params = 1:9, sim = 1:NumSims)

# Write a function to extract the initial dispersal values
DispFunc <- function(i){
     # Sort out the parameter combination and simulation under consideration
     Param <- TraitIndices$params[i]
     AllSims <- list.files(paste("~/ShiftingSlopes/MainSim/Params", Param, sep = ""))
     SimID <- AllSims[TraitIndices$sim[i]]
     
     # Load the corresponding initial population matrix
     InFile <- paste("~/ShiftingSlopes/MainSim/Params", Param, "/", 
                     SimID, "/InitialPopMat.csv", sep = "")
     SysCommand1 <- paste("gunzip ", InFile, ".gz", sep = "")
     SysCommand2 <- paste("gzip ", InFile, sep = "")
     system(SysCommand1)
     SimData <- read.csv(InFile)
     system(SysCommand2)
     
     # Calculate the dispersal phenotypes
     DispCols <- grep(pattern = "disp", x = names(SimData))
     DispSums <- rowSums(SimData[, DispCols])
     ExpDists <- (dmax * eta * exp(rho * DispSums)) / (1 + exp(rho * DispSums))
     MuInit <- mean(ExpDists)
     
     # Now remove the PopMats and load in the summary statistics to calculate
     #    the change in mean dispersal phenotypes brought on by climate change
     rm(SimData)
     SpeedFolders <- c("Slow", "MainSim", "Fast")
     EndShift <- c(LengthShift, BurnIn + LengthShift, LengthShift)
     EvolResults <- matrix(NA, nrow = 3, ncol = 2)
     for(v in 1:3){
          InFile <- paste("~/ShiftingSlopes/", SpeedFolders[v], "/Params", Param,
                          "/", SimID, "/SummaryStats.csv", sep = "")
          SysCommand1 <- paste("gunzip ", InFile, ".gz", sep = "")
          SysCommand2 <- paste("gzip ", InFile, sep = "")
          system(SysCommand1)
          SimData <- read.csv(InFile)
          system(SysCommand2)
          # Calculate a generational correction for the x coordinates to convert them
          #    to relative values
          SimData$Relx <- rep(NA, nrow(SimData))
          for(j in 1:nrow(SimData)){
               CurGen <- ifelse(v == 2, SimData$gen[j], SimData$gen[j] + BurnIn)
               SimData$Relx[j] <- SimData$x[j] - BetaCoord[CurGen]
          }
          # Extract the data for the starting generation and form a data frame for
          #    the results
          if(v !=2 ){
               StartInFile <- paste("~/ShiftingSlopes/MainSim/Params", Param,
                                    "/", SimID, "/SummaryStats.csv", sep = "")
               SysCommand1 <- paste("gunzip ", StartInFile, ".gz", sep = "")
               SysCommand2 <- paste("gzip ", StartInFile, sep = "")
               system(SysCommand1)
               StartData <- read.csv(StartInFile)
               system(SysCommand2)
               StartData <- subset(StartData, (gen == BurnIn) & (abund >= Threshold))
               StartData$Relx <- StartData$x - BetaCoord[BurnIn]
          } else{
               StartData <- subset(SimData, (gen == BurnIn) & (abund >= Threshold))
          }
          DispEvol <- data.frame(x = StartData$Relx, y = StartData$y, InitDisp = StartData$muDisp,
                                 FinalDisp = rep(NA, nrow(StartData)))
          
          # Walk through the data to fill in the FinalDisp and Dist columns
          for(j in 1:nrow(DispEvol)){
               PatchData <- subset(SimData, (Relx == DispEvol$x[j]) & (y == DispEvol$y[j]) & (abund >= Threshold))
               if(max(PatchData$gen) > EndShift[v]){
                    PostGens <- subset(PatchData, gen >= EndShift[v])
                    FinalGen <- min(PostGens$gen)
               } else{
                    FinalGen <- max(PatchData$gen)
               }
               FinalData <- subset(PatchData, gen == FinalGen)
               if(dim(FinalData)[1] == 1){
                    DispEvol$FinalDisp[j] <- FinalData$muDisp
               }
          }
          # Calculate the average change in dispersal across the landscape
          DeltaDisp <- DispEvol$FinalDisp - DispEvol$InitDisp
          MuDelta <- mean(DeltaDisp, na.rm = TRUE)
          # Check to see if the simulation went extinct at this speed of climate change
          GenExtinct <- max(SimData$gen) + 1
          Ext <- ifelse(GenExtinct < (MaxGen[v]), 0, 1)
          # Put it all together
          EvolResults[v,] <- c(MuDelta, Ext)
     }
     
     Results <- list(MuInit, EvolResults)
     return(Results)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("TraitIndices", "dmax", "eta", "rho", "BurnIn", "Threshold",
                    "LengthShift", "MaxGen", "BetaCoord"))

# Run the simulations
SimVec <- 1:nrow(TraitIndices)
SimDisp <- clusterApply(cl, x = SimVec, fun = DispFunc)


save(SimDisp, TraitIndices, file = "DispIdea.rdata")



