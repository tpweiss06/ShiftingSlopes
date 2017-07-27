# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# To test the MatFill function, (1) Make sure the output only includes occupied
#    patches (and includes all of them), (2) Each patch as the correct number
#    of offspring (passed in as an argument), (3) Make sure x and y coordinates
#    aren't getting mixed up, (4) Test that selfing occurs in monoecious but not
#    dioecious scenarios

AllCols <- PopMatColNames(nFit = 1, nDisp = 1, monoecious = TRUE)
FitCols <- grep("^fit", AllCols)
DispCols <- grep("^disp", AllCols)
PopMat <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                     PopSize = 10, BetaInit = 0, width = 5, FitInit = 1,
                     FitDiv = 1, DispInit = 1, DispDiv = 1)
TraitMat <- CalcTraits(population = 1:10, PopMat = PopMat, FitColumns = FitCols,
                       DispColumns = DispCols)
FitVals <- RelFit(LocalSel = 1, tau = 100, beta = 0, traits = TraitMat, PopMat = PopMat, 
                  individuals = 1:nrow(PopMat), omega = 2)
OccPatches <- unique(PopMat[,c("x1", "y1")])
RealizedNtp1 <- c(1, 2, 3, 4, 5)
U <- c(0.5, 0.5)
names(U) <- c("fit", "disp")
Vm <- c(0.1, 0.1)
names(Vm) <- c("fit", "disp")

MatFill(RealizedNtp1 = RealizedNtp1, ColumnNames = AllCols, OccPatches = OccPatches, 
        SexRatio = 0.5, RelFits = FitVals, U = U, Vm = Vm, FitColumns = FitCols, 
        DispColumns = DispCols, PopMat = PopMat)



writeLines("Pass", con = "ModelTesting/MatFill.txt")










