# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# To test the CalcTraits function, (1) use the Initialize function to create
#    an input matrix and ensure CalcTraits returns ouput of the proper
#    dimensionality, (2) Plot histograms of all the trait values for all
#    individuals. Since they are summing across normally distributed random
#    variables, they should conform to a normal distribution with a mean of
#    2*NumLoci*TraitMean and a standard deviation of sqrt(2*NumLoci*TraitDiv^2).
#    For the dispersal trait, first take the log of the trat values, then
#    compare to this distribution since the dispersal sum is exponentiated.

NumFitLoci <- 12
NumDispLoci <- 7
TrialPopSize <- 1000
FitMean <- 2
FitDiv <- 3
DispMean <- 0.2
DispDiv <- 6

ColNames <- PopMatColNames(nFit = NumFitLoci, nDisp = NumDispLoci, monoecious = TRUE)
FitCols <- grep("^fit", ColNames)
DispCols <- grep("^disp", ColNames)
TrialPop <- Initialize(ColumnNames = ColNames, FitCols = FitCols, DispCols = DispCols, 
                       PopSize = TrialPopSize, BetaInit = 0, width = 3, 
                       SexRatio = 0.5, FitInit = FitMean, FitDiv = FitDiv, 
                       DispInit = DispMean, DispDiv = DispDiv)

TrialTraits <- CalcTraits(population = 1:nrow(TrialPop), PopMat = TrialPop, 
                          FitColumns = FitCols, DispColumns = DispCols)

if(nrow(TrialTraits) == TrialPopSize & ncol(TrialTraits) == 3){
     Check1 <- TRUE
} else{
     Check1 <- FALSE
}


hist(TrialTraits[,"fit"], freq = FALSE)
xseq <- seq(0, 100, length.out = 1000)
lines(x = xseq, dnorm(x = xseq, mean = 2*NumFitLoci*FitMean, sd = sqrt(2*NumFitLoci*FitDiv^2)))

hist(log(TrialTraits[,"disp"]), freq = FALSE)
xseq <- seq(-60, 80, length.out = 1000)
lines(x = xseq, dnorm(x = xseq, mean = 2*NumDispLoci*DispMean, sd = sqrt(2*NumDispLoci*DispDiv^2)))

Check2 <- TRUE

if(Check1 & Check2){
     writeLines("Pass", con = "ModelTesting/CalcTraits.txt")
} else{
     writeLines("Fail", con = "ModelTesting/CalcTraits.txt")
}






