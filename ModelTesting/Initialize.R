# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# To test the Initialize function: (1) First use PopMatColNames() to generate
#    some of the necessary input (2) Check that the returned matrix is of the
#    correct dimensions (3) All x0 locations should be equal to BetaInit and
#    y0 should be uniformally distributed between 1 and the width argument
#    (4) Play with the SexRatio argument to ensure that it works as desired
#    (5) Check that the allele values are correctly normally distributed.

TestColumnNames <- PopMatColNames(nFit = 2, nDisp = 2, monoecious = FALSE)
FitColumns <- grep("^fit", TestColumnNames)
DispColumns <- grep("^disp", TestColumnNames)

TestPop1 <- Initialize(ColumnNames = TestColumnNames, FitCols = FitColumns, 
                       DispCols = DispColumns, PopSize = 1000, BetaInit = 0, 
                       width = 10, SexRatio = 0.5, FitInit = 0, FitDiv = 1, 
                       DispInit = 0, DispDiv = 1)
# Check the dimensions
nrow(TestPop1) == 1000
ncol(TestPop1) == length(TestColumnNames)
Check1 <- TRUE

# Check the starting locations
unique(TestPop1[,"x1"])
unique(TestPop1[,"y1"])
# If the breaks are not set manually, then the individuals in y1 patches of 1 and
#    2 get lumped together leading to a less than uniform seeming distribution.
hist(TestPop1[,"y1"], breaks = c(0,1,2,3,4,5,6,7,8,9,10))
Check2 <- TRUE

# Check the sex ratio
mean(TestPop1[,"sex"])
TestPop2 <- Initialize(ColumnNames = TestColumnNames, FitCols = FitColumns, 
                       DispCols = DispColumns, PopSize = 1000, BetaInit = 0, 
                       width = 10, SexRatio = 0.75, FitInit = 0, FitDiv = 1, 
                       DispInit = 0, DispDiv = 1)
mean(TestPop2[,"sex"])
TestPop3 <- Initialize(ColumnNames = TestColumnNames, FitCols = FitColumns, 
                       DispCols = DispColumns, PopSize = 1000, BetaInit = 0, 
                       width = 10, SexRatio = 0.25, FitInit = 0, FitDiv = 1, 
                       DispInit = 0, DispDiv = 1)
mean(TestPop3[,"sex"])
Check3 <- TRUE


# Check the distribution of allele frequencies
TestPop4 <- Initialize(ColumnNames = TestColumnNames, FitCols = FitColumns, 
                       DispCols = DispColumns, PopSize = 5000, BetaInit = 0, 
                       width = 10, SexRatio = 0.25, FitInit = 0, FitDiv = 1, 
                       DispInit = 5, DispDiv = 1)
TestPop5 <- Initialize(ColumnNames = TestColumnNames, FitCols = FitColumns, 
                       DispCols = DispColumns, PopSize = 5000, BetaInit = 0, 
                       width = 10, SexRatio = 0.25, FitInit = 0, FitDiv = 3, 
                       DispInit = 5, DispDiv = 0.25)
xseq <- seq(-10, 15, length.out = 1000)
hist(TestPop4[,FitColumns], freq = FALSE)
lines(x = xseq, y = dnorm(x = xseq, mean = 0, sd = 1))

hist(TestPop4[,DispColumns], freq = FALSE)
lines(x = xseq, y = dnorm(x = xseq, mean = 5, sd = 1))

hist(TestPop5[,FitColumns], freq = FALSE)
lines(x = xseq, y = dnorm(x = xseq, mean = 0, sd = 3))

hist(TestPop5[,DispColumns], freq = FALSE)
lines(x = xseq, y = dnorm(x = xseq, mean = 5, sd = 0.25))

Check4 <- TRUE

if(Check1 & Check2 & Check3 & Check4){
     writeLines("Pass", con = "ModelTesting/Initialize.txt")
} else{
     writeLines("Fail", con = "ModelTesting/Initialize.txt")
}
