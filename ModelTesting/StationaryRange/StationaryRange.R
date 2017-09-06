# Script to test for normal population growth within a single patch in the 
#    absence of dispersal. 

setwd("~/Desktop/RangeShifts/ShiftingSlopes/ModelTesting/StationaryRange/")
source("~/Desktop/RangeShifts/ShiftingSlopes/SimFunctions.R")

# First set all parameters
BetaInit <- 0
gamma <- 0.025
tau <- 276.25
LocalSel <- 1
omega <- 3
U <- c(0.02, 0.02)
Vm <- c(0.0004, 0.0004)
nFit <- 5
nDisp <- 5
R0 <- 1
K0 <- 10
width <- 10
kern <- "exp"
EnvGradType <- "K"
monoecious <- TRUE
BurnIn <- 100
BurnOut <- 0
LengthShift <- 0
ClimSpeed <- 0
InitPopSize <- 50
FitInit <- 0.05
FitDiv <- 0.025
DispInit <- 0
DispDiv <- 1
PatchScale <- 10
NumRands <- 1000000
SexRatio <- 0.5

parameters <- list(BetaInit, gamma, tau, LocalSel, omega, U, Vm, nFit, nDisp, R0,
                   K0, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                   LengthShift, ClimSpeed, InitPopSize, FitInit, FitDiv, DispInit,
                   DispDiv, PatchScale, NumRands, SexRatio)
names(parameters) <- c("BetaInit", "gamma", "tau", "LocalSel", "omega", "U", "Vm", 
                       "nFit", "nDisp", "R0", "K0", "width", "kern", "EnvGradType", 
                       "monoecious", "BurnIn", "BurnOut", "LengthShift", "ClimSpeed", 
                       "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                       "PatchScale", "NumRands", "SexRatio")

FullSim(parameters = parameters, parallel = FALSE)


# Graph the resulting abundance patterns through time as a single, animated graph
SumStats <- read.csv("Sim1/SummaryStats.csv")
head(SumStats)
range(SumStats$x)
plot(SumStats$x, SumStats$abund)
quantile(SumStats$x, probs = c(0, 0.05, 0.10, 0.15, 0.25, 0.75, 0.85, 0.90, 0.95, 1))
QuantSumStats <- subset(SumStats, (x >= -83) & (x <= 91))
plot(QuantSumStats$x, QuantSumStats$abund)

# The population definitely moved to fill the range. Next steps are to put some
#    simulations on the server for longer and then to play with this data to make
#    nice, animated heat maps or some other good visualization...

