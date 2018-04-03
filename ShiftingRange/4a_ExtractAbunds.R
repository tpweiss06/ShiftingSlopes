# This script will extract the abundance information from the stationary range
#    simulations. Specifically, this script will save the mean abundance through
#    time at each sector (i.e. lattice  column with identical environmental 
#    conditions), the variance in abundance values at each patch
#    (i.e. among simulation variance), and the variance among
#    patches within the same sector (i.e. within simulation variance).

# Set the speed for the current script
SpeedIndex <- 1

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
     InFile <- paste("~/ShiftingSlopes/ShiftingRange", SpeedWords[SpeedIndex], 
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
     return(Abunds)
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
     AbundVals[AbundIndices$params[i], AbundIndices$sim[i],,,] <- SimAbunds[[i]]
}

# Make another function to pass to the cluster
AbundProcess <- function(p){
     ParamSectorMean <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     ParamAmongVar <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     ParamWithinVar <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     
     for(i in 1:RangeExtent){
          for(j in 1:NumGens){
               PatchMeans <- colMeans(AbundVals[p,,j,i,], na.rm = TRUE)
               ParamSectorMean[i,j] <- mean(PatchMeans, na.rm = TRUE)
               
               PatchVars <- rep(NA, width)
               for(k in 1:width){
                    PatchVars[k] <- var(AbundVals[p,,j,i,k], na.rm = TRUE)
               }
               ParamAmongVar[i,j] <- mean(PatchVars, na.rm = TRUE)
               
               SimVars <- rep(NA, 100)
               for(k in 1:100){
                    SimVars[k] <- var(AbundVals[p,k,j,i,], na.rm = TRUE)
               }
               ParamWithinVar[i,j] <- mean(SimVars, na.rm = TRUE)
          }
     }
     TempList <- list(SectorMean = ParamSectorMean, AmongVar = ParamAmongVar,
                      WithinVar = ParamWithinVar)
     return(TempList)
}

# Now create the summary objects to be saved from this array
SectorMean <- array(NA, dim = c(9, RangeExtent, NumGens))
AmongVar <- array(NA, dim = c(9, RangeExtent, NumGens))
WithinVar <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     SimSummary <- AbundProcess(p)
     SectorMean[p,,] <- SimSummary$SectorMean
     AmongVar[p,,] <- SimSummary$AmongVar
     WithinVar[p,,] <- SimSummary$WithinVar
}

# Finally save the output
OutFile <- paste(SpeedWords[SpeedIndex], "ShiftingAbundResults.rdata", sep = "")
save(SectorMean, AmongVar, WithinVar, file = OutFile)


