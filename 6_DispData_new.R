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
load("ExtSimIDs.rdata")

# Create the vectors to be passed to the function in the data extraction
SimVec <- 1:27
SpeedWords <- c("Slow", "MainSim", "Fast")
SpeedVec <- rep(SpeedWords, 9)
ParamVec <- rep(1:9, each = 3)

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
EndShift <- c(LengthShift, BurnIn + LengthShift, LengthShift)
BetaCoord <- matrix(NA, nrow = 3, ncol = BurnIn + LengthShift + BurnOut)
for(i in 1:3){
     BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
                                eta = eta, v = ClimSpeeds[i]) / eta
     BetaCoord[i,] <- c(rep(BetaInit, BurnIn), BetaShift, rep(BetaShift[LengthShift], BurnOut))
}

DispExtract <- function(i){
     # Set the current speed and parameter combination and isolate the relevant
     #    simulation IDs
     SpeedWord <- SpeedVec[i]
     ParamCombo <- ParamVec[i]
     CurSims <- subset(SimIDs, Params == ParamCombo)
     
     # Source a representative parameter file for useful values
     ParamFile <- paste("~/ShiftingSlopes/MainSim/Params", ParamCombo, "/", CurSims$SimID[1],
                        "/parameters.R", sep = "")
     source(ParamFile)
     
     # Loop through each simulation to get the values and store them in the vectors
     ExtantInit <- NULL
     ExtinctInit <- NULL
     ExtantEvolve <- NULL
     ExtinctEvolve <- NULL
     for(s in 1:nrow(CurSims)){
          InitPopPath <- paste("~/ShiftingSlopes/MainSim/Params", ParamCombo, "/", CurSims$SimID[s], "/InitialPopMat.csv", sep = "")
          EndPath <- paste("~/ShiftingSlopes/", SpeedWord, "/Params", ParamCombo, "/", CurSims$SimID[s], "/SummaryStats.csv", sep = "")
          InitPop <- read.csv(InitPopPath)
          EndData <- read.csv(EndPath)
          
          # Calculate the initial dispersal phenotypes and store them appropriately
          if(s == 1){ # Only do this once
               DispCols <- grep(pattern = "disp", x = names(InitPop))
          }
          DispSums <- rowSums(InitPop[, DispCols])
          ExpDists <- (dmax * eta * exp(rho * DispSums)) / (1 + exp(rho * DispSums))
          if(SpeedWord == "Slow"){
               if(CurSims$Slow[s] == 1){
                    ExtantInit <- c(ExtantInit, ExpDists)
               }else{
                    ExtinctInit <- c(ExtinctInit, ExpDists)
               }
          } else if(SpeedWord = "MainSim"){
               if(CurSims$Moderate[s] == 1){
                    ExtantInit <- c(ExtantInit, ExpDists)
               }else{
                    ExtinctInit <- c(ExtinctInit, ExpDists)
               }
          } else{
               if(CurSims$Fast[s] == 1){
                    ExtantInit <- c(ExtantInit, ExpDists)
               }else{
                    ExtinctInit <- c(ExtinctInit, ExpDists)
               }
          }
          
          # Now calculate the evolved changes in dispersal beginning with
          #    calculating relative x coordinates
          EndData$Relx <- rep(NA, nrow(EndData))
          for(j in 1:nrow(EndData)){
               CurGen <- ifelse(v == 2, EndData$gen[j], EndData$gen[j] + BurnIn)
               EndData$Relx[j] <- EndData$x[j] - BetaCoord[CurGen]
          }
          # Isolate the start data we need
          if(SpeedWord == "MainSim"){
               StartData <- subset(EndData, (gen == BurnIn) & (abund >= Threshold))
          } else{
               InitDataPath <- paste("~/ShiftingSlopes/MainSim/Params", ParamCombo, "/", CurSims$SimID[s], "/SummaryStats.csv", sep = "")
               InitData <- read.csv(InitDataPath)
               StartData <- subset(InitData, (gen == BurnIn) & (abund >= Threshold))
               StartData$Relx <- StartData$x - BetaCoord[BurnIn]
          }
          DispEvol <- data.frame(x = StartData$Relx, y = StartData$y, InitDisp = StartData$muDisp,
                                 FinalDisp = rep(NA, nrow(StartData)))
          
          # Walk through the data to fill in the FinalDisp and Dist columns
          for(j in 1:nrow(DispEvol)){
               PatchData <- subset(EndData, (Relx == DispEvol$x[j]) & (y == DispEvol$y[j]) & (abund >= Threshold))
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
          # Now calculate the vector of changes in average dispersal phenotypes
          DeltaDisp <- DispEvol$FinalDisp - DispEvol$InitDisp
          # Finally, store it in the appropriate vector
          if(SpeedWord == "Slow"){
               if(CurSims$Slow[s] == 1){
                    ExtantEvolve <- c(ExtantEvolve, DeltaDisp)
               }else{
                    ExtinctEvolve <- c(ExtinctEvolve, DeltaDisp)
               }
          } else if(SpeedWord = "MainSim"){
               if(CurSims$Moderate[s] == 1){
                    ExtantEvolve <- c(ExtantEvolve, DeltaDisp)
               }else{
                    ExtinctEvolve <- c(ExtinctEvolve, DeltaDisp)
               }
          } else{
               if(CurSims$Fast[s] == 1){
                    ExtantEvolve <- c(ExtantEvolve, DeltaDisp)
               }else{
                    ExtinctEvolve <- c(ExtinctEvolve, DeltaDisp)
               }
          }
     }
     Results <- list(ExtantInit = ExtantInit, ExtinctInit = ExtinctInit,
                     ExtantEvolve = ExtantEvolve, ExtinctEvolve = ExtinctEvolve)
     return(Results)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("dmax", "eta", "rho", "BurnIn", "Threshold", "SpeedVec",
                    "LengthShift", "MaxGen", "BetaCoord", "ParamVec", "SimIDs"))

# Run the simulations
SimDisp <- clusterApply(cl, x = SimVec, fun = DispExtract)

# Process the results into appropriate objects for the later graphing script
DispInit <- vector(mode = "list", length = 9)
DispEvolve <- vector(mode = "list", length = 9)
for(p in 1:9){
     DispInit[[p]] <- vector(mode = "list", length = 3)
     DispEvolve[[p]] <- vector(mode = "list", length = 3)
     for(v in 1:3){
          CurIndex <- which((SpeedVec == SpeedWords[v]) & (ParamVec == p))
          CurResults <- SimFits[[CurIndex]]
          DispInit[[p]][[v]] <- list(Extant = CurResults$ExtantInit, 
                                     Extinct = CurResults$ExtinctInit)
          DispEvolve[[p]][[v]] <- list(Extant = CurResults$ExtantEvolve, 
                                     Extinct = CurResults$ExtinctEvolve)
     }
}

# Finally, save the output
save(DispInit, file = "~/ShiftingSlopes/DispInitNew.rdata")
save(DispEvolve, file = "~/ShiftingSlopes/DispEvolveNew.rdata")

