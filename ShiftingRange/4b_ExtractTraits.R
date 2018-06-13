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
     Success[TraitIndices$params[i], TraitIndices$sim[i]] <- !SimTraits[[i]]$Ext
}

# Make another function to pass to the cluster
TraitProcess <- function(p, FocalSims){
     MeanPhenotypeFit <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     MeanPhenotypeDisp <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     MeanPatchPhenCVFit <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     MeanPatchPhenCVDisp <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     MeanPatchGenSigmaFit <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     MeanPatchGenSigmaDisp <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     WithinPhenCVFit <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     WithinPhenCVDisp <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     AmongPhenCVFit <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     AmongPhenCVDisp <- matrix(NA, nrow = RangeExtent, ncol = NumGens)

     for(i in 1:RangeExtent){
          for(j in 1:NumGens){
               MeanPhenotypeFit[i,j] <- mean(FitVals[p,FocalSims,j,i,,1], na.rm = TRUE)
               MeanPhenotypeDisp[i,j] <- mean(DispVals[p,FocalSims,j,i,,1], na.rm = TRUE)
               
               MeanPatchPhenCVFit[i,j] <- mean(FitVals[p,FocalSims,j,i,,2] / FitVals[p,FocalSims,j,i,,1],
                                               na.rm = TRUE)
               MeanPatchPhenCVDisp[i,j] <- mean(DispVals[p,FocalSims,j,i,,2] / DispVals[p,FocalSims,j,i,,1],
                                               na.rm = TRUE)
               MeanPatchGenSigmaFit[i,j] <- mean(FitVals[p,FocalSims,j,i,,3], na.rm = TRUE)
               MeanPatchGenSigmaDisp[i,j] <- mean(DispVals[p,FocalSims,j,i,,3], na.rm = TRUE)
               
               AmongPhenCVFit[i,j] <- sd(FitVals[p,FocalSims,j,i,,1], na.rm = TRUE) /
                                        MeanPhenotypeFit[i,j]
               AmongPhenCVDisp <- sd(DispVals[p,FocalSims,j,i,,1], na.rm = TRUE) /
                                        MeanPhenotypeDisp[i,j]
               if(length(FocalSims) > 1){
                    FitCVs <- rep(NA, length(FocalSims))
                    DispCVs <- rep(NA, length(FocalSims))
                    for(k in 1:length(FocalSims)){
                         FitCVs[k] <- sd(FitVals[p,FocalSims[k],j,i,,1], na.rm = TRUE) / 
                              mean(FitVals[p,FocalSims[k],j,i,,1], na.rm = TRUE)
                         DispCVs[k] <- sd(DispVals[p,FocalSims[k],j,i,,1], na.rm = TRUE) / 
                              mean(DispVals[p,FocalSims[k],j,i,,1], na.rm = TRUE)
                    }
                    WithinPhenCVFit[i,j] <- mean(FitVars, na.rm = TRUE)
                    WithinPhenCVDisp[i,j] <- mean(DispVars, na.rm = TRUE)
               } else{
                    WithinPhenCVFit[i,j] <- ParamAmongVarFit[i,j]
                    WithinPhenCVDisp[i,j] <- ParamAmongVarDisp[i,j]
               }
          }
     }
     TempList <- list(MeanPhenotypeFit = MeanPhenotypeFit,
                      MeanPhenotypeDisp = MeanPhenotypeDisp,
                      MeanPatchPhenCVFit = MeanPatchPhenCVFit,
                      MeanPatchPhenCVDisp = MeanPatchPhenCVDisp,
                      MeanPatchGenSigmaFit = MeanPatchGenSigmaFit,
                      MeanPatchGenSigmaDisp = MeanPatchGenSigmaDisp,
                      WithinPhenCVFit = WithinPhenCVFit,
                      WithinPhenCVDisp = WithinPhenCVDisp,
                      AmongPhenCVFit = AmongPhenCVFit,
                      AmongPhenCVDisp = AmongPhenCVDisp)
     return(TempList)
}

# Now process the results for those simulations that successfully tracked climate
#    change
SuccessMeanPhenotypeFit <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessMeanPhenotypeDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessMeanPatchPhenCVFit <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessMeanPatchPhenCVDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessMeanPatchGenSigmaFit <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessMeanPatchGenSigmaDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessWithinPhenCVFit <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessWithinPhenCVDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessAmongPhenCVFit <- array(NA, dim = c(9, RangeExtent, NumGens))
SuccessAmongPhenCVDisp <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     ParamSims <- Success[p,]
     SuccessSims <- which(ParamSims == TRUE)
     if(length(SuccessSims) > 0){
          SimSummary <- TraitProcess(p, FocalSims = SuccessSims)
          SuccessMeanPhenotypeFit[p,,] <- SimSummary$MeanPhenotypeFit
          SuccessMeanPhenotypeDisp[p,,] <- SimSummary$MeanPhenotypeDisp
          SuccessMeanPatchPhenCVFit[p,,] <- SimSummary$MeanPatchPhenCVFit
          SuccessMeanPatchPhenCVDisp[p,,] <- SimSummary$MeanPatchPhenCVDisp
          SuccessMeanPatchGenSigmaFit[p,,] <- SimSummary$MeanPatchGenSigmaFit
          SuccessMeanPatchGenSigmaDisp[p,,] <- SimSummary$MeanPatchGenSigmaDisp
          SuccessWithinPhenCVFit[p,,] <- SimSummary$WithinPhenCVFit
          SuccessWithinPhenCVDisp[p,,] <- SimSummary$WithinPhenCVDisp
          SuccessAmongPhenCVFit[p,,] <- SimSummary$AmongPhenCVFit
          SuccessAmongPhenCVDisp[p,,] <- SimSummary$AmongPhenCVDisp
     }
}

# Now the simulations that failed to track climate change
FailureMeanPhenotypeFit <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureMeanPhenotypeDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureMeanPatchPhenCVFit <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureMeanPatchPhenCVDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureMeanPatchGenSigmaFit <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureMeanPatchGenSigmaDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureWithinPhenCVFit <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureWithinPhenCVDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureAmongPhenCVFit <- array(NA, dim = c(9, RangeExtent, NumGens))
FailureAmongPhenCVDisp <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     ParamSims <- Success[p,]
     FailureSims <- which(ParamSims == FALSE)
     if(length(FailureSims) > 0){
          SimSummary <- TraitProcess(p, FocalSims = FailureSims)
          FailureMeanPhenotypeFit[p,,] <- SimSummary$MeanPhenotypeFit
          FailureMeanPhenotypeDisp[p,,] <- SimSummary$MeanPhenotypeDisp
          FailureMeanPatchPhenCVFit[p,,] <- SimSummary$MeanPatchPhenCVFit
          FailureMeanPatchPhenCVDisp[p,,] <- SimSummary$MeanPatchPhenCVDisp
          FailureMeanPatchGenSigmaFit[p,,] <- SimSummary$MeanPatchGenSigmaFit
          FailureMeanPatchGenSigmaDisp[p,,] <- SimSummary$MeanPatchGenSigmaDisp
          FailureWithinPhenCVFit[p,,] <- SimSummary$WithinPhenCVFit
          FailureWithinPhenCVDisp[p,,] <- SimSummary$WithinPhenCVDisp
          FailureAmongPhenCVFit[p,,] <- SimSummary$AmongPhenCVFit
          FailureAmongPhenCVDisp[p,,] <- SimSummary$AmongPhenCVDisp
     }
}

# Now all of them
TotalMeanPhenotypeFit <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalMeanPhenotypeDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalMeanPatchPhenCVFit <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalMeanPatchPhenCVDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalMeanPatchGenSigmaFit <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalMeanPatchGenSigmaDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalWithinPhenCVFit <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalWithinPhenCVDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalAmongPhenCVFit <- array(NA, dim = c(9, RangeExtent, NumGens))
TotalAmongPhenCVDisp <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     ParamSims <- Success[p,]
     TotalSims <- 1:length(ParamSims)
     if(length(TotalSims) > 0){
          SimSummary <- TraitProcess(p, FocalSims = TotalSims)
          TotalMeanPhenotypeFit[p,,] <- SimSummary$MeanPhenotypeFit
          TotalMeanPhenotypeDisp[p,,] <- SimSummary$MeanPhenotypeDisp
          TotalMeanPatchPhenCVFit[p,,] <- SimSummary$MeanPatchPhenCVFit
          TotalMeanPatchPhenCVDisp[p,,] <- SimSummary$MeanPatchPhenCVDisp
          TotalMeanPatchGenSigmaFit[p,,] <- SimSummary$MeanPatchGenSigmaFit
          TotalMeanPatchGenSigmaDisp[p,,] <- SimSummary$MeanPatchGenSigmaDisp
          TotalWithinPhenCVFit[p,,] <- SimSummary$WithinPhenCVFit
          TotalWithinPhenCVDisp[p,,] <- SimSummary$WithinPhenCVDisp
          TotalAmongPhenCVFit[p,,] <- SimSummary$AmongPhenCVFit
          TotalAmongPhenCVDisp[p,,] <- SimSummary$AmongPhenCVDisp
     }
}

# Finally save the output
OutFile <- paste(SpeedWords[SpeedIndex], "ShiftingTraitResults.rdata", sep = "")
MeanPhenotypeFit <- list(Success = SuccessMeanPhenotypeFit, Failure = FailureMeanPhenotypeFit, 
                  Total = TotalMeanPhenotypeFit)
MeanPhenotypeDisp <- list(Success = SuccessMeanPhenotypeDisp, Failure = FailureMeanPhenotypeDisp, 
                          Total = TotalMeanPhenotypeDisp)
MeanPatchPhenCVFit <- list(Success = SuccessMeanPatchPhenCVFit, Failure = FailureMeanPatchPhenCVFit, 
                           Total = TotalMeanPatchPhenCVFit)
MeanPatchPhenCVDisp <- list(Success = SuccessPatchPhenCVDisp, Failure = FailurePatchPhenCVDisp, 
                            Total = TotalPatchPhenCVDisp)
MeanPatchGenSigmaFit <- list(Success = SuccessMeanPatchGenSigmaFit, Failure = FailureMeanPatchGenSigmaFit, 
                             Total = TotalMeanPatchGenSigmaFit)
MeanPatchGenSigmaDisp <- list(Success = SuccessMeanPatchGenSigmaDisp, Failure = FailureMeanPatchGenSigmaDisp, 
                              Total = TotalMeanPatchGenSigmaDisp)
WithinPhenCVFit <- list(Success = SuccessWithinPhenCVFit, Failure = FailureWithinPhenCVFit, 
                        Total = TotalWithinPhenCVFit)
WithinPhenCVDisp <- list(Success = SuccessWithinPhenCVDisp, Failure = FailureWithinPhenCVDisp, 
                         Total = TotalWithinPhenCVDisp)
AmongPhenCVFit <- list(Success = SuccessAmongPhenCVFit, Failure = FailureAmongPhenCVFit, 
                       Total = TotalAmongPhenCVFit)
AmongPhenCVDisp <- list(Success = SuccessAmongPhenCVDisp, Failure = FailureAmongPhenCVDisp, 
                        Total = TotalAmongPhenCVDisp)

# Save all the output
save(MeanPhenotypeFit, MeanPhenotypeDisp, MeanPatchPhenCVFit, MeanPatchPhenCVDisp, 
     MeanPatchGenSigmaFit, MeanPatchGenSigmaDisp, WithinPhenCVFit, WithinPhenCVDisp,
     AmongPhenCVFit, AmongPhenCVDisp, Success, file = OutFile)

