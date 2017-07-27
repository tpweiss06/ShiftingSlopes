# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# To test the RelFit function, (1) create a large population with all x locations
#    at beta and the fitness trait mean set at 0.5, which should always be the 
#    optimal trait value at x = beta. Use RelFit to calculate the relative
#    fitness values of the results and plot a histogram. Depending on the setting
#    of the breaks, there should be a peak at 1 or just below one as all values
#    will be high, but many of them will be offset from 0.5 in either direction,
#    leading to lower fitness values than 1. (2) Repeat the above process with 
#    larger initial variance in the trait values to ensure the histogram adjusts
#    accordingly. (3) Repeat the same process with the original trait variance 
#    but a smaller omega, the results should be similar to 2. (4) Alter the 
#    population matrix to have x values from accross the range and repeat check
#    the fitness values with the same LocalSel parameter. Then do it again with 
#    LocalSel set to 0. The first scenario should lead to a peak with symmetric
#    declines in fitness across the range and the second should lead to relatively
#    constant fitness across the range.

# (1) Basic test

AllCols <- PopMatColNames(nFit = 1, nDisp = 4, monoecious = TRUE)
FitCols <- grep("^fit", AllCols)
DispCols <- grep("^disp", AllCols)
PopMat <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                         PopSize = 1000, BetaInit = 0, width = 5, FitInit = 0.25,
                         FitDiv = 0.25, DispInit = 0, DispDiv = 0.2)
TraitMat <- CalcTraits(population = 1:1000, PopMat = PopMat, FitColumns = FitCols,
                       DispColumns = DispCols)

FitVals <- RelFit(LocalSel = 1, tau = 100, beta = 0, traits = TraitMat, PopMat = PopMat, 
                  individuals = 1:nrow(PopMat), omega = 2)
hist(FitVals)
#hist(TraitMat[,"fit"])

Check1 <- TRUE

# (2) Check that a larger initial variance leads to lower fitness overall
AllCols <- PopMatColNames(nFit = 1, nDisp = 4, monoecious = TRUE)
FitCols <- grep("^fit", AllCols)
DispCols <- grep("^disp", AllCols)
PopMat <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                     PopSize = 1000, BetaInit = 0, width = 5, FitInit = 0.25,
                     FitDiv = 1, DispInit = 0, DispDiv = 0.2)
TraitMat <- CalcTraits(population = 1:1000, PopMat = PopMat, FitColumns = FitCols,
                       DispColumns = DispCols)

FitVals <- RelFit(LocalSel = 1, tau = 100, beta = 0, traits = TraitMat, PopMat = PopMat, 
                  individuals = 1:nrow(PopMat), omega = 2)
hist(FitVals)

Check2 <- TRUE

# (3) Check that a lower omega leads to similar results
AllCols <- PopMatColNames(nFit = 1, nDisp = 4, monoecious = TRUE)
FitCols <- grep("^fit", AllCols)
DispCols <- grep("^disp", AllCols)
PopMat <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                     PopSize = 1000, BetaInit = 0, width = 5, FitInit = 0.25,
                     FitDiv = 0.25, DispInit = 0, DispDiv = 0.2)
TraitMat <- CalcTraits(population = 1:1000, PopMat = PopMat, FitColumns = FitCols,
                       DispColumns = DispCols)

FitVals <- RelFit(LocalSel = 1, tau = 100, beta = 0, traits = TraitMat, PopMat = PopMat, 
                  individuals = 1:nrow(PopMat), omega = 0.5)
hist(FitVals)

Check3 <- TRUE

# (4) Checkt that local selection is working
AllCols <- PopMatColNames(nFit = 1, nDisp = 4, monoecious = TRUE)
FitCols <- grep("^fit", AllCols)
DispCols <- grep("^disp", AllCols)
PopMat <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                     PopSize = 1000, BetaInit = 0, width = 5, FitInit = 0.25,
                     FitDiv = 0, DispInit = 0, DispDiv = 0.2)
PopMat[,"x1"] <- seq(-500, 500, length.out = 1000)
TraitMat <- CalcTraits(population = 1:1000, PopMat = PopMat, FitColumns = FitCols,
                       DispColumns = DispCols)

FitVals1 <- RelFit(LocalSel = 5, tau = 100, beta = 0, traits = TraitMat, PopMat = PopMat, 
                  individuals = 1:nrow(PopMat), omega = 0.25)
FitVals2 <- RelFit(LocalSel = 0, tau = 100, beta = 0, traits = TraitMat, PopMat = PopMat, 
                   individuals = 1:nrow(PopMat), omega = 0.25)
plot(x = PopMat[,"x1"], y = FitVals1, col = "darkred", xlab = "location", 
     ylab = "Relative fitness", main = "")
points(x = PopMat[,"x1"], y = FitVals2, col = "blue")

Check4 <- TRUE

if(Check1 & Check2 & Check3 & Check4){
     writeLines("Pass", con = "ModelTesting/RelFits.txt")
}else{
     writeLines("Fail", con = "ModelTesting/RelFits.txt")
}


