# This script will run simulations on the MSI compute cluster to characterize
#    population structure in a stable range for all of the 9 range parameter
#    combinations to be explored.

# Set the number of processors and number of simulations to be run
nProc <- 24*18
NumSims <- 200

# Set the working directory and load necessary data and libraries
source("~/ShiftingSlopes/SimFunctions.R")
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")
library(parallel)
library(Rmpi)

# Now load in the parameters and compile all 9 combinations into a single list
CurSpeed <- 0.5
AllParams <- vector(mode = "list", length = 9)
for(i in 1:9){
     CurDirectory <- paste("~/ShiftingSlopes/MainSim/Params", i, sep = "")
     ExampleSim <- list.files(CurDirectory)[1]
     ParamFile <- paste(CurDirectory, "/", ExampleSim, "/parameters.R", sep = "")
     source(ParamFile)
     
     v <- CurSpeed
     BurnIn <- 0
     AllParams[[i]] <- list(BetaInit, gamma, tau, lambda, omega, U, Vm, Lf, Ld, Rmax,
                         Kmax, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                         LengthShift, v, InitPopSize, FitInit, FitDiv, DispInit,
                         DispDiv, eta, NumRands, z, dmax, rho)
     names(AllParams[[i]]) <- c("BetaInit", "gamma", "tau", "lambda", "omega", "U", "Vm", 
                              "Lf", "Ld", "Rmax", "Kmax", "width", "kern", "EnvGradType", 
                              "monoecious", "BurnIn", "BurnOut", "LengthShift", "v", 
                              "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                              "eta", "NumRands", "z", "dmax", "rho")
}

# Write a function to be passed to various nodes
SimFunc <- function(i){
     # Load in the appropriate initial population matrix
     p <- ParamVec[i]
     InFile <- paste("~/ShiftingSlopes/MainSim/Params", p, "/", 
                     StatSimIDVec[i], "/InitialPopMat.csv", sep = "")
     InputMat <- read.csv(InFile)
     InputMat$x1 <- InputMat$x0
     InputMat$y1 <- InputMat$y0
     InputMat <- as.matrix(InputMat)
     setwd(paste("~/ShiftingSlopes/Slow/Params", p, "/", sep = ""))
     FullSim(parameters = AllParams[[p]], parallel = TRUE, PopInit = InputMat, SimID = StatSimIDVec[i])
     return(i)
}

# Create a vector of parameter index values for the parallel computation
ParamVec <- rep(1:9, each = NumSims)
StatSimIDVec <- NULL
for(p in 1:9){
     FilePath <- paste("~/ShiftingSlopes/MainSim/Params", p ,sep = "")
     SimFiles <- list.files(FilePath)
     StatSimIDVec <- c(StatSimIDVec, SimFiles)
}
SimVec <- 1:(NumSims*9)

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("AllParams", "SimVec", "ParamVec", "StatSimIDVec") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/ShiftingSlopes/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


