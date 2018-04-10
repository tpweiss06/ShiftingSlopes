# This script will extract the trait information from the stationary range
#    simulations. For both fitness and dispersal the mean value at each time
#    point and sector (i.e. lattice column with identical environmental 
#    conditions) for the phenotype, genetic standard deviation, and phenotypic
#    standard deviation. Additionally, the within and among simulation variance
#    will be calculated for the mean value of each trait.

# Set the speed for the current script
SpeedIndex <- 1

# Set the number of processors available for this scrip
nProc <- 2*24

# Set the working directory
setwd("~/ShiftingSlopes/ShiftingRange/")
source("~/ShiftingSlopes/SimFunctions.R")
library(parallel)
library(Rmpi)
source("ShiftParams.R")
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")

# Create arrays to hold the trait values from each simulation. Dimensions 
#    of these arrays are: parameter combination, simulation number, time, x axis,
#    y axis, summary value.
#    x and y axis correspond to the center coordinates of each patch in the
#    discrete lattice
#    summary value ranges from 1 to 3, corresponding to the following:
#         1 -- MuTrait        (Mean phenotypic value)
#         2 -- SigTraitPhen   (Phenotypic standard deviation)
#         3 -- SigTraitGen    (Genetic standard deviation)
#    RangeExtent and ZeroPos below correspond to the mapping of unbounded real 
#    number x values corresponding to patch centers to array indices
NumGens <- 200
width <- 10
RangeExtent <- 121
FitVals <- array(NA, dim = c(9, 100, NumGens, RangeExtent, width, 3))
DispVals <- array(NA, dim = c(9, 100, NumGens, RangeExtent, width, 3))
Success <- matrix(NA, nrow = 9, ncol = 100)

BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
                           eta = RangeParams$eta[1], v = SpeedNums[SpeedIndex])
BetaShift <- BetaShift / RangeParams$eta[1]
BetaCoord <- c(rep(BetaInit, BurnIn), BetaShift, rep(BetaShift[LengthShift], BurnOut))
ZeroPos <- 61

# Create a data frame to index each simulation block in the above arrays
TraitIndices <- expand.grid(params = 1:9, sim = 1:100)

# Write a function to be passed to various nodes
TraitExtract <- function(i){
     # Sort out the parameter combination and simulation under consideration
     Param <- TraitIndices$params[i]
     AllSims <- list.files(paste("~/ShiftingSlopes/ShiftingRange/", SpeedWords[SpeedIndex], 
                                 "/Params", Param, "/", sep = ""))
     SimID <- AllSims[TraitIndices$sim[i]]
     
     # Load the corresponding summary statistics
     InFile <- paste("~/ShiftingSlopes/ShiftingRange/", SpeedWords[SpeedIndex], 
                     "/Params", Param, "/", SimID, "/SummaryStats.csv", sep = "")
     SimData <- read.csv(InFile)
     
     # Create an array matching the dimensions of the trait value arrays to 
     #    store the trait values from this particular simulation
     Fit <- array(NA, dim = c(NumGens, RangeExtent, width, 3))
     Disp <- array(NA, dim = c(NumGens, RangeExtent, width, 3))
     
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
                              Fit[g, xArrInd, k, 1] <- CurPatch$muFit
                              Fit[g, xArrInd, k, 2] <- CurPatch$sigmaFitPhen
                              Fit[g, xArrInd, k, 3] <- CurPatch$sigmaFitGen
                              Disp[g, xArrInd, k, 1] <- CurPatch$muDisp
                              Disp[g, xArrInd, k, 2] <- CurPatch$sigmaDispPhen
                              Disp[g, xArrInd, k, 3] <- CurPatch$sigmaDispGen
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
     results <- list(Fit = Fit, Disp = Disp, Ext = Ext)
     return(results)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("SpeedWords", "NumGens", "width", "RangeExtent", "ZeroPos", 
                    "TraitIndices", "BetaCoord", "SpeedIndex"))

# Run the simulations
SimVec <- 1:dim(TraitIndices)[1]
SimTraits <- clusterApply(cl, x = SimVec, fun = TraitExtract)

# Now populate the AbundVals array with the results
for(i in SimVec){
     FitVals[TraitIndices$params[i], TraitIndices$sim[i],,,,] <- SimTraits[[i]]$Fit
     DispVals[TraitIndices$params[i], TraitIndices$sim[i],,,,] <- SimTraits[[i]]$Disp
     Success[AbundIndices$params[i], AbundIndices$sim[i]] <- !SimAbunds[[i]]$Ext
}

# Make another function to pass to the cluster
TraitProcess <- function(p, FocalSims){
     ParamSectorFit <- array(NA, dim = c(RangeExtent, NumGens, 3))
     ParamSectorDisp <- array(NA, dim = c(RangeExtent, NumGens, 3))
     ParamAmongVarFit <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     ParamAmongVarDisp <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     ParamWithinVarFit <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     ParamWithinVarDisp <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     
     for(i in 1:RangeExtent){
          for(j in 1:NumGens){
               for(v in 1:3){
                    FitMeans <- colMeans(FitVals[p,FocalSims,j,i,,v], na.rm = TRUE)
                    ParamSectorFit[i,j,v] <- mean(FitMeans, na.rm = TRUE)
                    
                    DispMeans <- colMeans(DispVals[p,FocalSims,j,i,,v], na.rm = TRUE)
                    ParamSectorDisp[i,j,v] <- mean(DispMeans, na.rm = TRUE)
               }
               
               FitVars <- rep(NA, width)
               DispVars <- rep(NA, width)
               for(k in 1:width){
                    FitVars[k] <- var(FitVals[p,FocalSims,j,i,k,1], na.rm = TRUE)
                    DispVars[k] <- var(DispVals[p,FocalSims,j,i,k,1], na.rm = TRUE)
               }
               ParamAmongVarFit[i,j] <- mean(FitVars, na.rm = TRUE)
               ParamAmongVarDisp[i,j] <- mean(DispVars, na.rm = TRUE)
               
               FitVars <- rep(NA, length(FocalSims))
               DispVars <- rep(NA, length(FocalSims))
               for(k in 1:length(FocalSims)){
                    FitVars[k] <- var(FitVals[p,FocalSims[k],j,i,,1], na.rm = TRUE)
                    DispVars[k] <- var(DispVals[p,FocalSims[k],j,i,,1], na.rm = TRUE)
               }
               ParamWithinVarFit[i,j] <- mean(FitVars, na.rm = TRUE)
               ParamWithinVarDisp[i,j] <- mean(DispVars, na.rm = TRUE)
          }
     }
     TempList <- list(SectorFit = ParamSectorFit, SectorDisp = ParamSectorDisp, AmongVarDisp = ParamAmongVarDisp,
                      AmongVarFit = ParamAmongVarFit, WithinVarDisp = ParamWithinVarDisp,
                      WithinVarFit = ParamWithinVarFit)
     return(TempList)
}

# Now process the results for those simulations that successfully tracked climate
#    change
SuccessSectorFit <- array(NA, dim = c(9, RangeExtent, NumGens, 3))
SuccessSectorDisp <- array(NA, dim = c(9, RangeExtent, NumGens, 3))
SuccessAmongVarFit <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessAmongVarDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessWithinVarFit <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessWithinVarDisp <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     ParamSims <- Success[p,]
     SimSummary <- TraitProcess(p, FocalSims = which(ParamSims == TRUE))
     SuccessSectorFit[p,,,] <- SimSummary$SectorFit
     SuccessSectorDisp[p,,,] <- SimSummary$SectorDisp
     SuccessAmongVarFit[p,,] <- SimSummary$AmongVarFit
     SuccessAmongVarDisp[p,,] <- SimSummary$AmongVarDisp
     SuccessWithinVarFit[p,,] <- SimSummary$WithinVarFit
     SuccessWithinVarDisp[p,,] <- SimSummary$WithinVarDisp
}

# Now the simulations that failed to track climate change
FailureSectorFit <- array(NA, dim = c(9, RangeExtent, NumGens, 3))
FailureSectorDisp <- array(NA, dim = c(9, RangeExtent, NumGens, 3))
FailureAmongVarFit <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureAmongVarDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureWithinVarFit <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureWithinVarDisp <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     ParamSims <- Success[p,]
     SimSummary <- TraitProcess(p, FocalSims = which(ParamSims == FALSE))
     FailureSectorFit[p,,,] <- SimSummary$SectorFit
     FailureSectorDisp[p,,,] <- SimSummary$SectorDisp
     FailureAmongVarFit[p,,] <- SimSummary$AmongVarFit
     FailureAmongVarDisp[p,,] <- SimSummary$AmongVarDisp
     FailureWithinVarFit[p,,] <- SimSummary$WithinVarFit
     FailureWithinVarDisp[p,,] <- SimSummary$WithinVarDisp
}


# Now all of them
TotalSectorFit <- array(NA, dim = c(9, RangeExtent, NumGens, 3))
TotalSectorDisp <- array(NA, dim = c(9, RangeExtent, NumGens, 3))
TotalAmongVarFit <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalAmongVarDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalWithinVarFit <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalWithinVarDisp <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     ParamSims <- Success[p,]
     SimSummary <- TraitProcess(p, FocalSims = 1:length(ParamSims))
     TotalSectorFit[p,,,] <- SimSummary$SectorFit
     TotalSectorDisp[p,,,] <- SimSummary$SectorDisp
     TotalAmongVarFit[p,,] <- SimSummary$AmongVarFit
     TotalAmongVarDisp[p,,] <- SimSummary$AmongVarDisp
     TotalWithinVarFit[p,,] <- SimSummary$WithinVarFit
     TotalWithinVarDisp[p,,] <- SimSummary$WithinVarDisp
}

# Finally save the output
OutFile <- paste(SpeedWords[SpeedIndex], "ShiftingTraitResults.rdata", sep = "")
SectorFit <- list(Success = SuccessSectorFit, Failure = FailureSectorFit, 
                  Total = TotalSectorFit)
AmongVarFit <- list(Success = SuccessAmongVarFit, Failure = FailureAmongVarFit, 
                    Total = TotalAmongVarFit)
WithinVarFit <- list(Success = SuccessWithinVarFit, Failure = FailureWithinVarFit, 
                      Total = TotalWithinVarFit)
SectorDisp <- list(Success = SuccessSectorDisp, Failure = FailureSectorDisp, 
                   Total = TotalSectorDisp)
AmongVarDisp <- list(Success = SuccessAmongVarDisp, Failure = FailureAmongVarDisp, 
                     Total = TotalAmongVarDisp)
WithinVarDisp <- list(Success = SuccessWithinVarDisp, Failure = FailureWithinVarDisp, 
                      Total = TotalWithinVarDisp)

# Save all the output
save(SectorFit, AmongVarFit, WithinVarFit, SectorDisp, AmongVarDisp, WithinVarDisp,
     file = OutFile)


