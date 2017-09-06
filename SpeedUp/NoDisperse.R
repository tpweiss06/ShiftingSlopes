# Script to test for normal population growth within a single patch in the 
#    absence of dispersal. 

setwd("~/Desktop/RangeShifts/ShiftingSlopes/ModelTesting/NoDisperse/")
source("~/Desktop/RangeShifts/ShiftingSlopes/SimulationFunctions.R")

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
BurnIn <- 5
BurnOut <- 0
LengthShift <- 0
ClimSpeed <- 0
InitPopSize <- 50
FitInit <- 0.05
FitDiv <- 0.025
DispInit <- -100
DispDiv <- 0
PatchScale <- 10

parameters <- list(BetaInit, gamma, tau, LocalSel, omega, U, Vm, nFit, nDisp, R0,
                   K0, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                   LengthShift, ClimSpeed, InitPopSize, FitInit, FitDiv, DispInit,
                   DispDiv, PatchScale)
names(parameters) <- c("BetaInit", "gamma", "tau", "LocalSel", "omega", "U", "Vm", 
                       "nFit", "nDisp", "R0", "K0", "width", "kern", "EnvGradType", 
                       "monoecious", "BurnIn", "BurnOut", "LengthShift", "ClimSpeed", 
                       "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                       "PatchScale")

FullSim(parameters = parameters, parallel = FALSE)


# Graph the resulting abundance patterns
quartz(width = 7, height = 5)
par(mfrow = c(2, 5), mar = c(0, 0, 0, 0), oma = c(5, 4, 4, 2) + 0.1)
load("Sim1/SummaryStats.Rdata")
OccPatches <- rep(NA, width)
for(i in 1:width){
     PopData <- PatchAbund[,i,]
     OccPatches[i] <- which(PopData[])
}




