# This script will extract the trait information from the stationary range
#    simulations. For both fitness and dispersal the mean value at each time
#    point and sector (i.e. lattice column with identical environmental 
#    conditions) for the phenotype, genetic standard deviation, and phenotypic
#    standard deviation. Additionally, the within and among simulation variance
#    will be calculated for the mean value of each trait.

# Set the number of processors available for this scrip
nProc <- 24*18

# Set the working directory
setwd("~/ShiftingSlopes/StationaryRange/")
library(Rmpi)
library(parallel)

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
MaxGen <- 2000
width <- 10
RangeExtent <- 121
ZeroPos <- 61
FitVals <- array(NA, dim = c(9, 100, NumGens, RangeExtent, width, 3))
DispVals <- array(NA, dim = c(9, 100, NumGens, RangeExtent, width, 3))

# Create a data frame to index each simulation block in the above arrays
TraitIndices <- expand.grid(params = 1:9, sim = 1:100)

# Write a function to be passed to various nodes
TraitExtract <- function(i){
     # Sort out the parameter combination and simulation under consideration
     Param <- TraitIndices$params[i]
     AllSims <- list.files(paste("Params", Param, "/", sep = ""))
     SimID <- AllSims[TraitIndices$sim[i]]
     
     # Load the corresponding summary statistics
     InFile <- paste("Params", Param, "/", SimID, "/SummaryStats.csv", sep = "")
     SimData <- read.csv(InFile)
     
     # Create an array matching the dimensions of the trait value arrays to 
     #    store the trait values from this particular simulation
     Fit <- array(NA, dim = c(NumGens, RangeExtent, width, 3))
     Disp <- array(NA, dim = c(NumGens, RangeExtent, width, 3))
     
     for(g in (MaxGen - NumGens + 1):MaxGen){
          CurGen <- subset(SimData, (gen == g) & (abund > 0))
          if(dim(CurGen)[1] > 0){
               xRange <- range(CurGen$x)
               xSeq <- seq(xRange[1], xRange[2], by = 1)
               for(j in xSeq){
                    CurCol <- subset(CurGen, x == j)
                    xArrInd <- ZeroPos + j
                    if(xArrInd <= RangeExtent){
                         for(k in 1:width){
                              CurPatch <- subset(CurCol, y == k)
                              if(dim(CurPatch)[1] == 1){
                                   Fit[g - (MaxGen - NumGens), xArrInd, k, 1] <- CurPatch$muFit
                                   Fit[g - (MaxGen - NumGens), xArrInd, k, 2] <- CurPatch$sigmaFitPhen
                                   Fit[g - (MaxGen - NumGens), xArrInd, k, 3] <- CurPatch$sigmaFitGen
                                   Disp[g - (MaxGen - NumGens), xArrInd, k, 1] <- CurPatch$muDisp
                                   Disp[g - (MaxGen - NumGens), xArrInd, k, 2] <- CurPatch$sigmaDispPhen
                                   Disp[g - (MaxGen - NumGens), xArrInd, k, 3] <- CurPatch$sigmaDispGen
                              }
                         }
                    }     
               }
          } 
     }
     results <- list(Fit = Fit, Disp = Disp)
     return(results)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("NumGens", "MaxGen", "width", "RangeExtent", "ZeroPos", "TraitIndices"))
temp <- clusterEvalQ(cl, setwd("~/ShiftingSlopes/StationaryRange/"))

# Run the simulations
SimVec <- 1:dim(TraitIndices)[1]
SimTraits <- clusterApply(cl, x = SimVec, fun = TraitExtract)

# Now populate the AbundVals array with the results
for(i in SimVec){
     FitVals[TraitIndices$params[i], TraitIndices$sim[i],,,,] <- SimTraits[[i]]$Fit
     DispVals[TraitIndices$params[i], TraitIndices$sim[i],,,,] <- SimTraits[[i]]$Disp
}

# Make another function to pass to the cluster
TraitProcess <- function(p){
     ParamSectorFit <- array(NA, dim = c(RangeExtent, NumGens, 3))
     ParamSectorDisp <- array(NA, dim = c(RangeExtent, NumGens, 3))
     ParamAmongVarFit <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     ParamWithinVarDisp <- matrix(NA, nrow = RangeExtent, ncol = NumGens)
     
     for(i in 1:RangeExtent){
          for(j in 1:NumGens){
               for(v in 1:3){
                    FitMeans <- colMeans(FitVals[p,,j,i,,v], na.rm = TRUE)
                    ParamSectorFit[i,j,v] <- mean(FitMeans, na.rm = TRUE)
                    
                    DispMeans <- colMeans(DispVals[p,,j,i,,v], na.rm = TRUE)
                    ParamSectorDisp[i,j,v] <- mean(DispMeans, na.rm = TRUE)
               }
               
               FitVars <- rep(NA, width)
               DispVars <- rep(NA, width)
               for(k in 1:width){
                    FitVars[k] <- var(FitVals[p,,j,i,k,1], na.rm = TRUE)
                    DispVars[k] <- var(DispVals[p,,j,i,k,1], na.rm = TRUE)
               }
               ParamAmongVarFit[i,j] <- mean(FitVars, na.rm = TRUE)
               ParamAmongVarDisp[i,j] <- mean(DispVars, na.rm = TRUE)
               
               FitVars <- rep(NA, 100)
               DispVars <- rep(NA, 100)
               for(k in 1:100){
                    FitVars[k] <- var(FitVals[p,k,j,i,,1], na.rm = TRUE)
                    DispVars[k] <- var(DispVals[p,k,j,i,,1], na.rm = TRUE)
               }
               ParamWithinVarFit[i,j] <- mean(FitVars, na.rm = TRUE)
               ParamWithinVarDisp[i,j] <- mean(DispVars, na.rm = TRUE)
          }
     }
     TempList <- list(SectorFit = ParamSectorFit, SectorDisp = ParamSectorDisp,
                      AmongVarFit = ParamAmongVarFit, WithinVarDisp = ParamWithinVarDisp)
     return(TempList)
}

clusterExport(cl, c("FitVals", "DispVals"))

SimSummary <- clusterApply(cl, x = 1:9, fun = TraitProcess)

# Now create the summary objects to be saved from these arrays
SectorFit <- array(NA, dim = c(9, RangeExtent, NumGens, 3))
SectorDisp <- array(NA, dim = c(9, RangeExtent, NumGens, 3))
AmongVarFit <- array(NA, dim = c(9, RangeExtent, NumGens))
AmongVarDisp <- array(NA, dim = c(9, RangeExtent, NumGens))
WithinVarFit <- array(NA, dim = c(9, RangeExtent, NumGens))
WithinVarDisp <- array(NA, dim = c(9, RangeExtent, NumGens))

for(p in 1:9){
     SectorFit[p,,,] <- SimSummary[[p]]$SectorFit
     SectorDisp[p,,,] <- SimSummary[[p]]$SectorDisp
     AmongVarFit[p,,] <- SimSummary[[p]]$AmongVarFit
     AmongVarDisp[p,,] <- SimSummary[[p]]$AmongVarDisp
     WithinVarFit[p,,] <- SimSummary[[p]]$WithinVarFit
     WithinVarDisp[p,,] <- SimSummary[[p]]$WithinVarDisp
}

# Finally save the output
save(SectorFit, AmongVarFit, WithinVarFit, SectorDisp, AmongVarDisp, WithinVarDisp,
     file = "StationaryTraitResults.rdata")


