# This script will run simulations on the MSI compute cluster to characterize
#    population structure in a stable range for one of the 9 range parameter
#    combinations that will be explored.

# Set the number of processors and number of simulations to be run
nProc <- 24
NumSims <- round((nProc - 1) * 4, digits = -2)


# Set the working directory and load necessary data and libraries
setwd("~/ShiftingSlopes/StationaryRange/")
source("~/ShiftingSlopes/SimFunctions.R")
RangeParams <- read.csv("RangeParameters.csv")
library(parallel)
library(Rmpi)

# Set a parameter index to keep track of which combination of range parameters
#    is currently being explored.
ParamIndex <- 2
ParamFolder <- "~/ShiftingSlopes/StationaryRange/Params2/"

# Now set all parameters
BetaInit <- 0
omega <- 3
U <- c(0.02, 0.02)
Vm <- c(0.0004, 0.0004)
nFit <- 5
nDisp <- 5
R0 <- 1
K0 <- 1
width <- 10
kern <- "exp"
EnvGradType <- "K"
monoecious <- TRUE
BurnIn <- 1000
BurnOut <- 0
LengthShift <- 0
ClimSpeed <- 0
InitPopSize <- 50
FitInit <- 0.05
FitDiv <- 0.025
DispInit <- 0
DispDiv <- 1
PatchScale <- 100
NumRands <- 1000000
SexRatio <- 0.5
gamma <- RangeParams$gamma[ParamIndex]
LocalSel <- RangeParams$LocalSel[ParamIndex]
tau <- RangeParams$tau[ParamIndex]

parameters <- list(BetaInit, gamma, tau, LocalSel, omega, U, Vm, nFit, nDisp, R0,
                   K0, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                   LengthShift, ClimSpeed, InitPopSize, FitInit, FitDiv, DispInit,
                   DispDiv, PatchScale, NumRands, SexRatio)
names(parameters) <- c("BetaInit", "gamma", "tau", "LocalSel", "omega", "U", "Vm", 
                       "nFit", "nDisp", "R0", "K0", "width", "kern", "EnvGradType", 
                       "monoecious", "BurnIn", "BurnOut", "LengthShift", "ClimSpeed", 
                       "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                       "PatchScale", "NumRands", "SexRatio")

# Write a function to be passed to various nodes
SimFunc <- function(i){
     FullSim(parameters = parameters, parallel = TRUE)
     return(i)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Check the cluster for a response from each node
# clusterCall( cl, function() Sys.info()[c("nodename","machine")])

# Export the necessary objects to each node
clusterExport(cl, c("parameters", "ParamFolder", "NumSims") )

# Change the working directory of the worker nodes
clusterEvalQ(cl, setwd(ParamFolder) )
clusterEvalQ(cl, source("~/ShiftingSlopes/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = 1:NumSims, fun = SimFunc)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


