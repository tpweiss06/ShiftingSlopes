# Script to test for normal population growth within a single patch in the 
#    absence of dispersal. 

# First set all parameters
BetaInit <- 0
gamma <- 0.025
tau <- 276.25
LocalSel <- 1
omega <- 3
U <- c(0.02, 0)
Vm <- c(0.0004, 0)
names(U) <- names(Vm) <- c("fit", "disp")
nFit <- 5
nDisp <- 5
R0 <- 1
K0 <- 10
width <- 10
kern <- "exp"
EnvGradType <- "K"
monoecious <- TRUE
BurnIn <- 20
BurnOut <- 0
LengthShift <- 0
ClimSpeed <- 0
InitPopSize <- 50
FitInit <- 0.05
FitDiv <- 0.025
DispInit <- -100
DispDiv <- 0
PatchScale <- 10
NumRands <- 1000000
SexRatio <- 0.5

parameters1 <- list(BetaInit, gamma, tau, LocalSel, omega, U, Vm, nFit, nDisp, R0,
                   K0, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                   LengthShift, ClimSpeed, InitPopSize, FitInit, FitDiv, DispInit,
                   DispDiv, PatchScale, NumRands, SexRatio)
names(parameters1) <- c("BetaInit", "gamma", "tau", "LocalSel", "omega", "U", "Vm", 
                       "nFit", "nDisp", "R0", "K0", "width", "kern", "EnvGradType", 
                       "monoecious", "BurnIn", "BurnOut", "LengthShift", "ClimSpeed", 
                       "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                       "PatchScale", "NumRands", "SexRatio")

parameters2 <- list(BetaInit, gamma, tau, LocalSel, omega, U, Vm, nFit, nDisp, R0,
                   K0, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                   LengthShift, ClimSpeed, InitPopSize, FitInit, FitDiv, DispInit,
                   DispDiv, PatchScale)
names(parameters2) <- c("BetaInit", "gamma", "tau", "LocalSel", "omega", "U", "Vm", 
                       "nFit", "nDisp", "R0", "K0", "width", "kern", "EnvGradType", 
                       "monoecious", "BurnIn", "BurnOut", "LengthShift", "ClimSpeed", 
                       "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                       "PatchScale")

library(microbenchmark)
SpedUpFunction <- function(parameters){
     source("~/Desktop/RangeShifts/ShiftingSlopes/SpeedUp/SimFunctionsSpeedUp.R")
     FullSim(parameters = parameters, parallel = FALSE)
}
OldFunction <- function(parameters){
     source("~/Desktop/RangeShifts/ShiftingSlopes/SimulationFunctions.R")
     FullSim(parameters = parameters, parallel = FALSE)
}

microbenchmark(SpedUpFunction(parameters = parameters1), OldFunction(parameters = parameters2))
#expr          min       lq        mean      median         uq        max       neval     cld
#SpedUp        675.7097  731.8522  774.9079  756.3489       810.5408  952.7566  100       a 
#OldFunction   1345.1347 1495.7278 1574.2803 1564.1136      1629.1451 1965.6205 100       b


