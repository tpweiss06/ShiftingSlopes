# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# To test the Inheritence function, I will use small, manually created input
#    population matrices to test specific cases. (1) Test cases with low and
#    high mutation rates and variance to ensure that the function behaves as
#    expected. (2) Set mutation to 0 and make a small input matrix with all
#    unique allele values to ensure that segregation is working correctly.

AllCols <- PopMatColNames(nFit = 1, nDisp = 4, monoecious = TRUE)
FitCols <- grep("^fit", AllCols)
DispCols <- grep("^disp", AllCols)
PopMat <- matrix(NA, nrow = 6, ncol = length(AllCols))
colnames(PopMat) <- AllCols
PopMat[,5:6] <- 1:12
PopMat[,7:14] <- 1:48

parents <- matrix(1:6, nrow = 3, ncol = 2)

Inheritence(Cols = FitCols, parents = parents, PopMat = PopMat, U = c(0.5, 0.5),
            Vm = c(1, 1))
Inheritence(Cols = FitCols, parents = parents, PopMat = PopMat, U = c(0.5, 0.5),
            Vm = c(0, 0))
Inheritence(Cols = FitCols, parents = parents, PopMat = PopMat, U = c(0, 0),
            Vm = c(1, 1))


Inheritence(Cols = DispCols, parents = parents, PopMat = PopMat, U = c(0.5, 0.5),
            Vm = c(1, 1))
Inheritence(Cols = DispCols, parents = parents, PopMat = PopMat, U = c(0.5, 0.5),
            Vm = c(0, 0))
Inheritence(Cols = DispCols, parents = parents, PopMat = PopMat, U = c(0, 0),
            Vm = c(1, 1))

writeLines("Pass", con = "ModelTesting/Inheritence.txt")
