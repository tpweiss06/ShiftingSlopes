# This script will generate the starting population matrices for the range
#    shifting simulations for all range parameter combinations

# Define the range parameter values that will be used and the number of processors
ParamVec <- 1:9
nProc <- 10

library(parallel)
library(Rmpi)

# Load the extracted trait and abundance data and set the zero position used in
#    the extraction scripts
load("StationaryTraitResults.rdata")
load("StationaryAbundResults.rdata")
ZeroPos <- 61

# Create a function to generate an expected population matrix for each parameter
#    combination
GenInputMat <- function(p){
     # Source the corresponding parameter file
     source(paste("~/ShiftingSlopes/ShiftingRange/ParamFiles/Params", p, ".R", sep = ""))
     
     # Create the population column indices
     PopIndices <- PopMatColNames(Lf = Lf, Ld = Ld, monoecious = monoecious)
     
     # Next generate the values to use in generating the matrix from the last
     #    generation of stationarity
     LastGen <- dim(SectorMean)[3]
     MeanSectorAbunds <- round(SectorMean[p,,LastGen])
     MeanSectorFitAllele <- SectorFit[p,,LastGen,1] / (2*Lf)
     SigmaSectorFitAllele <- SectorFit[p,,LastGen,3]
     # For an explanation of the derivation of the MeanDispAllele approximation, 
     #    see model documentation
     MeanDist <- SectorDisp[p,,LastGen,1]
     MeanSectorDispAllele <- log(MeanDist / (dmax * eta - MeanDist)) / (rho * Ld)
     SigmaSectorDispAllele <- SectorDisp[p,,LastGen,3]

     # Now generate a vector of population abundances per patch
     TotalPatches <- length(MeanSectorAbunds)*width
     PatchAbunds <- rep(NA, TotalPatches)
     
     # When determining abundances, also keep track of the sector and width associated
     #    with each individual.
     SectorIndex <- NULL
     WidthIndex <- NULL
     
     for(i in 1:length(MeanSectorAbunds)){
          Indices <- ((i-1)*width + 1):(i*width)
          PopSizes <- rep(MeanSectorAbunds[i], width)
          PatchAbunds[Indices] <- PopSizes
          SectorIndex <- c(SectorIndex, rep(i, sum(PopSizes)))
          for(j in 1:width){
               WidthIndex <- c(WidthIndex, rep(j, PopSizes[j]))
          }
     }
     
     # Now initialize the matrix with an entry for each individual
     nCol <- ifelse(monoecious, max(PopIndices$DispCols), PopIndices$sex)
     InputMat <- matrix(NA, nrow = sum(MeanSectorAbunds), ncol = nCol)
     
     # Now step through the input matrix and fill in all the necessary information for
     #    each individual. 
     for(i in 1:nrow(InputMat)){
          InputMat[i, PopIndices$x1] <- SectorIndex[i] - ZeroPos
          InputMat[i, PopIndices$y1] <- WidthIndex[i]
          InputMat[i, PopIndices$FitCols] <- rnorm(n = Lf * 2, mean = MeanSectorFitAllele[SectorIndex[i]],
                                                   sd = SigmaSectorFitAllele[SectorIndex[i]])
          InputMat[i, PopIndices$DispCols] <- rnorm(n = Ld * 2, mean = MeanSectorDispAllele[SectorIndex[i]],
                                                    sd = SigmaSectorDispAllele[SectorIndex[i]])
          if(!monoecious){
               InputMat[i, PopIndices$sex] <- rbinom(n = 1, size = 1, prob = z)
          }
     }
     
     OutFile <- paste("~/ShiftingSlopes/ShiftingRange/InputMats/InputMat", p, ".rdata", sep = "")
     save(InputMat, file = OutFile)
     return(NULL)
}

# In the cluster I will need to source the simulation functions for each node
cl <- makeCluster(nProc - 1, type = "MPI")

# Load the simulation functions on all the nodes
temp <- clusterEvalQ(cl, source("~/ShiftingSlopes/SimFunctions.R") )

# Export the necessary data arrays
clusterExport(cl, c("SectorMean", "SectorFit", "SectorDisp", "ZeroPos"))

# Run the simulations
Sims <- clusterApply(cl, x = ParamVec, fun = GenInputMat)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.



