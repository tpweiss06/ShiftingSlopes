# This script will run a sanity check with a single parameter combination and
#    climate change with a speed of 0 to ensure there are no extinctions and if
#    there are that there are no intital condition constraints (I guess. This is a 
#    pretty stupid waste of time).

# Set the number of processors and number of simulations to be run
nProc <- 24*2
NumSims <- 50

# Set the working directory and load necessary data and libraries
setwd("~/ShiftingSlopes/ShiftingRange/")
source("~/ShiftingSlopes/SimFunctions.R")
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")
library(parallel)
library(Rmpi)

# Now make a parameter list for just a single set of parameters
CurParam <- 1
v <- 0
ParamFile <- paste("ParamFiles/Params", CurParam, ".R", sep = "")
source(ParamFile)
source("ShiftParams.R")
AllParams <- list(BetaInit, gamma, tau, lambda, omega, U, Vm, Lf, Ld, Rmax,
                  Kmax, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                  LengthShift, v, InitPopSize, FitInit, FitDiv, DispInit,
                  DispDiv, eta, NumRands, z, dmax, rho)
names(AllParams) <- c("BetaInit", "gamma", "tau", "lambda", "omega", "U", "Vm", 
                      "Lf", "Ld", "Rmax", "Kmax", "width", "kern", "EnvGradType", 
                      "monoecious", "BurnIn", "BurnOut", "LengthShift", "v", 
                      "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                      "eta", "NumRands", "z", "dmax", "rho")

# Write a function to be passed to various nodes
SimFunc <- function(i){
     InFile <- paste("~/ShiftingSlopes/StationaryRange/Params", CurParam, "/", 
                     StatSimIDVec[i], "/PopMat.csv", sep = "")
     SysCommand1 <- paste("gunzip ", InFile, ".gz", sep = "")
     SysCommand2 <- paste("gzip ", InFile, sep = "")
     system(SysCommand1)
     InputMat <- read.csv(InFile)
     system(SysCommand2)
     InputMat$x1 <- InputMat$x0
     InputMat$y1 <- InputMat$y0
     InputMat <- as.matrix(InputMat)
     NewDirectory <- "~/ShiftingSlopes/SanityCheck/"
     setwd(NewDirectory)
     FullSim(parameters = AllParams, parallel = TRUE, PopInit = InputMat)
     return(i)
}

# Create a vector of parameter index values for the parallel computation
StatSimVec <- list.files(paste("~/ShiftingSlopes/StationaryRange/Params", CurParam ,sep = ""))[1:NumSims]
SimVec <- 1:NumSims

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("AllParams", "SimVec", "CurParam", "StatSimIDVec") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/ShiftingSlopes/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


