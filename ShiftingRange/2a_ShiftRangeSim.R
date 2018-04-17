# This script will run simulations on the MSI compute cluster to characterize
#    population structure in a stable range for all of the 9 range parameter
#    combinations to be explored.

# Set the number of processors and number of simulations to be run
nProc <- 24*18
NumSims <- 100

# Set the working directory and load necessary data and libraries
setwd("~/ShiftingSlopes/ShiftingRange/")
source("~/ShiftingSlopes/SimFunctions.R")
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")
library(parallel)
library(Rmpi)

# Now load in the parameters and compile all 9 combinations into a single list
#    with sublists for the three speeds of climate change
AllParams <- vector(mode = "list", length = 9)
for(i in 1:9){
     ParamFile <- paste("ParamFiles/Params", i, ".R", sep = "")
     source(ParamFile)
     source("ShiftParams.R")
     
     AllParams[[i]] <- vector(mode = "list", length = 3)
     names(AllParams[[i]]) <- SpeedWords
     for(j in 1:3){
          v <- SpeedNums[j]
          AllParams[[i]][[j]] <- list(BetaInit, gamma, tau, lambda, omega, U, Vm, Lf, Ld, Rmax,
                                 Kmax, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                                 LengthShift, v, InitPopSize, FitInit, FitDiv, DispInit,
                                 DispDiv, eta, NumRands, z, dmax, rho)
          names(AllParams[[i]][[j]]) <- c("BetaInit", "gamma", "tau", "lambda", "omega", "U", "Vm", 
                                     "Lf", "Ld", "Rmax", "Kmax", "width", "kern", "EnvGradType", 
                                     "monoecious", "BurnIn", "BurnOut", "LengthShift", "v", 
                                     "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                                     "eta", "NumRands", "z", "dmax", "rho")
     }
}

# Write a function to be passed to various nodes
SimFunc <- function(index){
     j <- nchar(ParamVec[index])
     i <- round(ParamVec[index] / 10^(j-1), digits = 1)
     InFile <- paste("~/ShiftingSlopes/StationaryRange/Params", i, "/", 
                     StatSimIDVec[index], "/PopMat.csv", sep = "")
     InputMat <- read.csv(InFile)
     InputMat$x1 <- InputMat$x0
     InputMat$y1 <- InputMat$y0
     Params <- AllParams[[i]][[j]]
     NewDirectory <- paste("~/ShiftingSlopes/ShiftingRange/", SpeedWords[j], "/Params", i, "/", sep = "")
     setwd(NewDirectory)
     FullSim(parameters = Params, parallel = TRUE, PopInit = InputMat)
     return(i)
}

# Create a vector of parameter index values for the parallel computation
ParamVec <- rep(c(1:9, 1:9*10, 1:9*100), each = NumSims)
StatSimIDVec <- NULL
for(speed in 1:3){
     for(p in 1:9){
          FilePath <- paste("~/ShiftingSlopes/StationaryRange/Params", p ,sep = "")
          SimFiles <- list.files(FilePath)
          StatSimIDVec <- c(StatSimIDVec, SimFiles)
     }
}
SimVec <- 1:(NumSims*9*3)
     
# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("AllParams", "SimVec", "SpeedWords", "ParamVec", "StatSimIDVec") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/ShiftingSlopes/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


