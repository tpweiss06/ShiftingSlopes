# This script will run simulations on the MSI compute cluster to characterize
#    population structure in a stable range for one of the 9 range parameter
#    combinations that will be explored.

# Set the number of processors and number of simulations to be run
nProc <- 24
NumSims <- 100

# Set the working directory and load necessary data and libraries
setwd("~/ShiftingSlopes/ShiftingRange/")
source("~/ShiftingSlopes/SimFunctions.R")
library(parallel)
library(Rmpi)

# Set a parameter index to keep track of which combination of range parameters
#    is currently being explored and load the relevant parameters
ParamIndex <- 1
ExampleSim <- list.files(path = paste("~/ShiftingSlopes/StationaryRange/Params", ParamIndex, "/", sep = ""))[1]
ParamFile <- paste("~/ShiftingSlopes/StationaryRange/Params", ParamIndex, "/", 
                   ExampleSim, "/parameters.R", sep = "")
source(ParamFile)

# Adjust the time parameters to simulate climate change rather than stationarity
BurnIn <- 0
LengthShift <- 100
BurnOut <- 100
ClimSpeed <- PatchScale*2

# Now use those parameters to make the input list
parameters <- list(BetaInit, gamma, tau, LocalSel, omega, U, Vm, nFit, nDisp, R0,
                   K0, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                   LengthShift, ClimSpeed, InitPopSize, FitInit, FitDiv, DispInit,
                   DispDiv, PatchScale, NumRands, SexRatio)
names(parameters) <- c("BetaInit", "gamma", "tau", "LocalSel", "omega", "U", "Vm", 
                       "nFit", "nDisp", "R0", "K0", "width", "kern", "EnvGradType", 
                       "monoecious", "BurnIn", "BurnOut", "LengthShift", "ClimSpeed", 
                       "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                       "PatchScale", "NumRands", "SexRatio")

# Load the results from the stationary range simulations and recalculate the Beta
#    index
InputTraits <- paste("Params", ParamIndex, "ExtractedTraits.rdata", sep = "")
InputAbunds <- paste("Params", ParamIndex, "ExtractedAbunds.rdata", sep = "")
load(InputTraits)
load(InputAbunds)
BetaPos <- nrow(SectorMeans) %/% 2

# Use these to generate an input population matrix of the correct format
# First generate the number of columns needed
# Next, create the population column indices
PopIndices <- PopMatColNames(nFit = nFit, nDisp = nDisp, monoecious = monoecious)

# Next generate the values to use in generating the matrix by averaging the last
#    10 time points from the stationary range simulations to estimate equilibrium
#    values.
FinalTimeBlock <- (ncol(SectorMeans)-10):ncol(SectorMeans)
MeanSectorAbunds <- rowMeans(SectorMeans[,FinalTimeBlock])
SigmaSectorAbunds <- rowMeans(SectorSigma[,FinalTimeBlock])
MeanSectorMuFit <- rowMeans(SectorMuFit[,FinalTimeBlock], na.rm = TRUE)
MeanSectorSigFit <- rowMeans(SectorSigFit[,FinalTimeBlock], na.rm = TRUE)
MeanSectorMuDisp <- rowMeans(SectorMuDisp[,FinalTimeBlock], na.rm = TRUE)
MeanSectorSigDisp <- rowMeans(SectorSigDisp[,FinalTimeBlock], na.rm = TRUE)

# Now generate a vector of population abundances per patch
TotalPatches <- length(MeanSectorAbunds)*width
RealizedPatchAbunds <- rep(NA, TotalPatches)
# When determining abundances, also keep track of the sector and width associated
#    with each individual.
SectorIndex <- NULL
WidthIndex <- NULL
for(i in 1:length(MeanSectorAbunds)){
     Indices <- ((i-1)*width + 1):(i*width)
     PopSizes <- round(rnorm(n = width, mean = MeanSectorAbunds[i], sd = SigmaSectorAbunds[i]))
     PopSizes[PopSizes<0] <- 0
     RealizedPatchAbunds[Indices] <- PopSizes
     SectorIndex <- c(SectorIndex, rep(i, sum(PopSizes)))
     for(j in 1:width){
          WidthIndex <- c(WidthIndex, rep(j, PopSizes[j]))
     }
}

# Now initialize the matrix with an entry for each individual
InputMat <- matrix(NA, nrow = sum(RealizedPatchAbunds), ncol = max(PopIndices$DispCols))

# Now step through the input matrix and fill in all the necessary information for
#    each individual. 
for(i in 1:nrow(InputMat)){
     InputMat[i, PopIndices$x1] <- SectorIndex[i] - BetaPos
     InputMat[i, PopIndices$y1] <- WidthIndex[i]
     InputMat[i, PopIndices$FitCols] <- rnorm(n = nFit * 2, mean = MeanSectorMuFit[SectorIndex[i]],
                                              sd = MeanSectorSigFit[SectorIndex[i]])
     InputMat[i, PopIndices$DispCols] <- rnorm(n = nDisp * 2, mean = MeanSectorMuDisp[SectorIndex[i]],
                                               sd = MeanSectorSigDisp[SectorIndex[i]])
     if(!monoecious){
          InputMat[i, PopIndices$sex] <- rbinom(n = 1, size = 1, prob = 0.5)
     }
}

# Write a function to be passed to various nodes
SimFunc <- function(i){
     FullSim(parameters = parameters, parallel = TRUE, PopInit = InputPopMat)
     return(i)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("parameters", "ParamFolder", "NumSims", "InputPopMat") )

# Change the working directory of the worker nodes
ParamFolder <- paste("~/ShiftingSlopes/ShiftingRange/Params", ParamIndex, "/", sep = "")
clusterEvalQ(cl, setwd(ParamFolder) )
clusterEvalQ(cl, source("~/ShiftingSlopes/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = 1:NumSims, fun = SimFunc)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


