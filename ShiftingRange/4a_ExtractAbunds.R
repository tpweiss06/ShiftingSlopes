# This script will extract the abundance information from the stationary range
#    simulations. Specifically, this script will save the mean abundance through
#    time at each sector (i.e. lattice  column with identical environmental 
#    conditions), the variance in abundance values at each patch
#    (i.e. among simulation variance), and the variance among
#    patches within the same sector (i.e. within simulation variance).

# Set the speed for the current script
SpeedIndex <- 3

# Set the number of processors available for this script
nProc <- 2*24

# Set the working directory
setwd("~/ShiftingSlopes/ShiftingRange/")
source("~/ShiftingSlopes/SimFunctions.R")
library(parallel)
library(Rmpi)
source("ShiftParams.R")
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")

# Create an array to hold the abundance values from each simulation. Dimensions 
#    of this array are: parameter combination, simulation number, time, x axis,
#    y axis.
#    x and y axis correspond to the center coordinates of each patch in the
#    discrete lattice
#    RangeExtent and ZeroPos below correspond to the mapping of unbounded real 
#    number x values corresponding to patch centers to array indices
NumGens <- 200
width <- 10
RangeExtent <- 121
AbundVals <- array(NA, dim = c(9, 100, NumGens, RangeExtent, width))
Success <- matrix(NA, nrow = 9, ncol = 100)

BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
                           eta = RangeParams$eta[1], v = SpeedNums[SpeedIndex])
BetaShift <- BetaShift / RangeParams$eta[1]
BetaCoord <- c(rep(BetaInit, BurnIn), BetaShift, rep(BetaShift[LengthShift], BurnOut))
ZeroPos <- 61

# Create a data frame to index each simulation block in the above array
AbundIndices <- expand.grid(params = 1:9, sim = 1:100)

# Write a function to be passed to various nodes
AbundExtract <- function(i){
     # Sort out the parameter combination and simulation under consideration
     Param <- AbundIndices$params[i]
     AllSims <- list.files(paste("~/ShiftingSlopes/ShiftingRange/", SpeedWords[SpeedIndex], 
                                 "/Params", Param, "/", sep = ""))
     SimID <- AllSims[AbundIndices$sim[i]]
     
     # Load the corresponding summary statistics
     InFile <- paste("~/ShiftingSlopes/ShiftingRange/", SpeedWords[SpeedIndex], 
                     "/Params", Param, "/", SimID, "/SummaryStats.csv", sep = "")
     SimData <- read.csv(InFile)
     
     # Create an array matching the dimensions of AbundVals to store the 
     #    abundances from this particular simulation
     Abunds <- array(NA, dim = c(NumGens, RangeExtent, width))
     
     for(g in 1:NumGens){
          CurGen <- subset(SimData, (gen == g) & (abund > 0))
          if(dim(CurGen)[1] > 0){
               xRange <- range(CurGen$x)
               xSeq <- seq(xRange[1], xRange[2], by = 1)
               for(j in xSeq){
                    CurCol <- subset(CurGen, x == j)
                    xArrInd <- ZeroPos + (j - BetaCoord[g])
                    for(k in 1:width){
                         CurPatch <- subset(CurCol, y == k)
                         if(dim(CurPatch)[1] == 1){
                              Abunds[g, xArrInd, k] <- CurPatch$abund
                         }
                    }
               }
          } 
     }
     GenExtinct <- max(SimData$gen) + 1
     if(GenExtinct < (NumGens + 1)){
          Ext <- TRUE
     } else{
          Ext <- FALSE
     }
     return(list(Abunds = Abunds, Ext = Ext))
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("SpeedWords", "NumGens", "width", "RangeExtent", "ZeroPos", 
                    "AbundIndices", "BetaCoord", "SpeedIndex"))

# Run the simulations
SimVec <- 1:dim(AbundIndices)[1]
SimAbunds <- clusterApply(cl, x = SimVec, fun = AbundExtract)

# Now populate the AbundVals array with the results
for(i in SimVec){
     AbundVals[AbundIndices$params[i], AbundIndices$sim[i],,,] <- SimAbunds[[i]]$Abunds
     Success[AbundIndices$params[i], AbundIndices$sim[i]] <- !SimAbunds[[i]]$Ext
}

# Make another function to pass to the cluster
AbundProcess <- function(p, FocalSims){
     ParamSectorMean <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     ParamAmongVar <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     ParamWithinVar <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     
     for(i in 1:RangeExtent){
          for(j in 1:NumGens){
               if(length(FocalSims) > 1){
                    PatchMeans <- colMeans(AbundVals[p,FocalSims,j,i,], na.rm = TRUE)
                    ParamSectorMean[i,j] <- mean(PatchMeans, na.rm = TRUE)
               
                    PatchVars <- rep(NA, width)
                    for(k in 1:width){
                         PatchVars[k] <- var(AbundVals[p,FocalSims,j,i,k], na.rm = TRUE)
                    }
                    ParamAmongVar[i,j] <- mean(PatchVars, na.rm = TRUE)
               
                    SimVars <- rep(NA, length(FocalSims))
                    for(k in 1:length(FocalSims)){
                         SimVars[k] <- var(AbundVals[p,FocalSims[k],j,i,], na.rm = TRUE)
                    }
                    ParamWithinVar[i,j] <- mean(SimVars, na.rm = TRUE)
               } else{
                    ParamSectorMean[i,j] <- mean(AbundVals[p,FocalSims,j,i,], na.rm = TRUE)
                    ParamAmongVar[i,j] <- NA
                    ParamWithinVar[i,j] <- var(AbundVals[p,FocalSims,j,i,], na.rm = TRUE)
               }
          }
     }
     TempList <- list(SectorMean = ParamSectorMean, AmongVar = ParamAmongVar,
                      WithinVar = ParamWithinVar)
     return(TempList)
}

# Now process the results for those simulations that successfully tracked climate
#    change
SuccessSectorMean <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessAmongVar <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessWithinVar <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     ParamSims <- Success[p,]
     SuccessSims <- which(ParamSims == TRUE)
     if(length(SuccessSims > 0)){
          SimSummary <- AbundProcess(p, FocalSims = SuccessSims)
          SuccessSectorMean[p,,] <- SimSummary$SectorMean
          SuccessAmongVar[p,,] <- SimSummary$AmongVar
          SuccessWithinVar[p,,] <- SimSummary$WithinVar
     }
}

# Now those that failed to track climate change
FailureSectorMean <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureAmongVar <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureWithinVar <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     ParamSims <- Success[p,]
     SimSummary <- AbundProcess(p, FocalSims = which(ParamSims == FALSE))
     FailureSectorMean[p,,] <- SimSummary$SectorMean
     FailureAmongVar[p,,] <- SimSummary$AmongVar
     FailureWithinVar[p,,] <- SimSummary$WithinVar
}

# Now all of them
TotalSectorMean <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalAmongVar <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalWithinVar <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     ParamSims <- Success[p,]
     SimSummary <- AbundProcess(p, FocalSims = 1:length(ParamSims))
     TotalSectorMean[p,,] <- SimSummary$SectorMean
     TotalAmongVar[p,,] <- SimSummary$AmongVar
     TotalWithinVar[p,,] <- SimSummary$WithinVar
}

# Finally save the output
OutFile <- paste(SpeedWords[SpeedIndex], "ShiftingAbundResults.rdata", sep = "")
SectorMean <- list(Success = SuccessSectorMean, Failure = FailureSectorMean,
                   Total = TotalSectorMean)
AmongVar <- list(Success = SuccessAmongVar, Failure = FailureAmongVar,
                   Total = TotalAmongVar)
WithinVar <- list(Success = SuccessWithinVar, Failure = FailureWithinVar,
                   Total = TotalWithinVar)
save(SectorMean, AmongVar, WithinVar, file = OutFile)
