# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# To test the SaveParams() function of the model, create a list of parameter
#    values to save, run the function, and check that it creates an executable
#    file that will load in values identical to those created initially.

OriginalParams <- list(beta = 1, gamma = 2, tau = 3, LocalSel = 4, omega = 5,
                       U = c(6,7), Vm = c(8,9), nFit = 10, nDisp = 11, R0 = 12,
                       K0 = 13, width = 14, kern = 15, EnvGradType = 16,
                       monoecious = 17, BurnIn = 18, BurnOut = 19, LengthShift = 20,
                       ClimSpeed = 21, InitPopSize = 22, BetaInit = 23, FitInit = 24,
                       FitDiv = 25, DispInit = 26, DispDiv = 27)

SaveParams(parameters = OriginalParams, FilePath = "ModelTesting")
source("ModelTesting/parameters.R")

if( (OriginalParams$beta == beta) &
    (OriginalParams$gamma == gamma) &
    (OriginalParams$tau == tau) &
    (OriginalParams$LocalSel == LocalSel) &
    (OriginalParams$omega == omega) &
    (OriginalParams$U[1] == U[1]) &
    (OriginalParams$U[2] == U[2]) &
    (OriginalParams$Vm[1] == Vm[1]) &
    (OriginalParams$Vm[2] == Vm[2]) &
    (OriginalParams$nFit == nFit) &
    (OriginalParams$nDisp == nDisp) &
    (OriginalParams$R0 == R0) &
    (OriginalParams$K0 == K0) &
    (OriginalParams$width == width) &
    (OriginalParams$kern == kern) &
    (OriginalParams$EnvGradType == EnvGradType) &
    (OriginalParams$monoecious == monoecious) &
    (OriginalParams$BurnIn == BurnIn) &
    (OriginalParams$BurnOut == BurnOut) &
    (OriginalParams$LengthShift == LengthShift) &
    (OriginalParams$ClimSpeed == ClimSpeed) &
    (OriginalParams$InitPopSize == InitPopSize) &
    (OriginalParams$BetaInit == BetaInit) &
    (OriginalParams$FitInit == FitInit) &
    (OriginalParams$FitDiv == FitDiv) &
    (OriginalParams$DispInit == DispInit) &
    (OriginalParams$DispDiv == DispDiv)){
     # Write out a pass file
     writeLines("Pass", con = "ModelTesting/SaveParams.txt")
} else{
     # Write out a fail file
     writeLines("Fail", con = "ModelTesting/SaveParams.txt")
}

system("rm ModelTesting/parameters.R")