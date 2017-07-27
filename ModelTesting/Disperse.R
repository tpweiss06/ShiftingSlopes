# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# To test the dispersal function, first use the previously vetted PopMatColNames,
#    Initialize, and CalcTraits functions to generate the necessary input. Then,
#    (1) set the landscape width to 5 with a large population size and use the
#    dispersal function to generate post dispersal locations according to
#    all three dispersal kernels. Examine them with a heat map and confirm
#    that the kernels have some of the expected differences in patterns.
#    (2) Repeat the above process with larger values of the dispersal trait in
#    the Intialize function and confirm that it leads to greater dispersal in
#    all three kernels.
#    (3) For a single kernel, keep all other parameters the same, but set the
#    width to 5, 10, and 20 with population sizes that will result to equivalent
#    average abundances across the width of the landscape (e.g. 1,000, 10,000, and
#    20,000).

library(plot3D)

AllCols <- PopMatColNames(nFit = 4, nDisp = 4, monoecious = TRUE)
FitCols <- grep("^fit", AllCols)
DispCols <- grep("^disp", AllCols)
PreDispMat <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                         PopSize = 1000, BetaInit = 0, width = 5, FitInit = 0,
                         FitDiv = 1, DispInit = 0, DispDiv = 0.2)
PreDispMat[,"x0"] <- PreDispMat[,"x1"]
PreDispMat[,"y0"] <- PreDispMat[,"y1"]

TraitMat <- CalcTraits(population = 1:1000, PopMat = PreDispMat, FitColumns = FitCols,
                       DispColumns = DispCols)

NormDisp <- Disperse(PopMat = PreDispMat, traits = TraitMat, width = 5, kern = "norm",
                     PatchScale = 1)
ExpDisp <- Disperse(PopMat = PreDispMat, traits = TraitMat, width = 5, kern = "exp",
                    PatchScale = 1)
StudDisp <- Disperse(PopMat = PreDispMat, traits = TraitMat, width = 5, kern = "stud_t",
                     PatchScale = 1)

NormMin <- min(NormDisp[,"x1"])
NormMax <- max(NormDisp[,"x1"])
NormHeatMap <- matrix(NA, nrow = 5, ncol = NormMax - NormMin + 1)
for(x in 0:(NormMax - NormMin)){
     for(y in 1:5){
          NormHeatMap[y,x+1] <- length(which((NormDisp[,"x1"] == (NormMin + x)) &
                                             NormDisp[,"y1"] == y))
     }
}
image2D(t(NormHeatMap), x = NormMin:NormMax, y = 1:5, main = "Normal Kernel")

ExpMin <- min(ExpDisp[,"x1"])
ExpMax <- max(ExpDisp[,"x1"])
ExpHeatMap <- matrix(NA, nrow = 5, ncol = ExpMax - ExpMin + 1)
for(x in 0:(ExpMax - ExpMin)){
     for(y in 1:5){
          ExpHeatMap[y,x+1] <- length(which((ExpDisp[,"x1"] == (ExpMin + x)) &
                                                 ExpDisp[,"y1"] == y))
     }
}
image2D(t(ExpHeatMap), x = ExpMin:ExpMax, y = 1:5, main = "Exponential Kernel")

StudMin <- min(StudDisp[,"x1"])
StudMax <- max(StudDisp[,"x1"])
StudHeatMap <- matrix(NA, nrow = 5, ncol = StudMax - StudMin + 1)
for(x in 0:(StudMax - StudMin)){
     for(y in 1:5){
          StudHeatMap[y,x+1] <- length(which((StudDisp[,"x1"] == (StudMin + x)) &
                                                  StudDisp[,"y1"] == y))
     }
}
image2D(t(StudHeatMap), x = StudMin:StudMax, y = 1:5, main = "Students't Kernel")

Check1 <- TRUE
# By adjusting values in the above code, I confirmed (2) as well
Check2 <- TRUE

# Adjust widths
PreDispMat1 <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                         PopSize = 1000, BetaInit = 0, width = 5, FitInit = 0,
                         FitDiv = 1, DispInit = 0, DispDiv = 0.2)
PreDispMat1[,"x0"] <- PreDispMat1[,"x1"]
PreDispMat1[,"y0"] <- PreDispMat1[,"y1"]

PreDispMat2 <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                          PopSize = 2000, BetaInit = 0, width = 10, FitInit = 0,
                          FitDiv = 1, DispInit = 0, DispDiv = 0.2)
PreDispMat2[,"x0"] <- PreDispMat2[,"x1"]
PreDispMat2[,"y0"] <- PreDispMat2[,"y1"]

PreDispMat3 <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                          PopSize = 4000, BetaInit = 0, width = 20, FitInit = 0,
                          FitDiv = 1, DispInit = 0, DispDiv = 0.2)
PreDispMat3[,"x0"] <- PreDispMat3[,"x1"]
PreDispMat3[,"y0"] <- PreDispMat3[,"y1"]

TraitMat1 <- CalcTraits(population = 1:1000, PopMat = PreDispMat1, FitColumns = FitCols,
                       DispColumns = DispCols)
TraitMat2 <- CalcTraits(population = 1:2000, PopMat = PreDispMat2, FitColumns = FitCols,
                        DispColumns = DispCols)
TraitMat3 <- CalcTraits(population = 1:4000, PopMat = PreDispMat3, FitColumns = FitCols,
                        DispColumns = DispCols)

ExpDisp1 <- Disperse(PopMat = PreDispMat1, traits = TraitMat1, width = 5, kern = "exp",
                     PatchScale = 1)
ExpDisp2 <- Disperse(PopMat = PreDispMat2, traits = TraitMat2, width = 10, kern = "exp",
                     PatchScale = 1)
ExpDisp3 <- Disperse(PopMat = PreDispMat3, traits = TraitMat3, width = 20, kern = "exp",
                     PatchScale = 1)

ExpMin1 <- min(ExpDisp1[,"x1"])
ExpMax1 <- max(ExpDisp1[,"x1"])
ExpHeatMap1 <- matrix(NA, nrow = 5, ncol = ExpMax1 - ExpMin1 + 1)
for(x in 0:(ExpMax1 - ExpMin1)){
     for(y in 1:5){
          ExpHeatMap1[y,x+1] <- length(which((ExpDisp1[,"x1"] == (ExpMin1 + x)) &
                                                 ExpDisp1[,"y1"] == y))
     }
}
image2D(t(ExpHeatMap1), x = ExpMin1:ExpMax1, y = 1:5, main = "Width = 5")


ExpMin2 <- min(ExpDisp2[,"x1"])
ExpMax2 <- max(ExpDisp2[,"x1"])
ExpHeatMap2 <- matrix(NA, nrow = 5, ncol = ExpMax2 - ExpMin2 + 1)
for(x in 0:(ExpMax2 - ExpMin2)){
     for(y in 1:5){
          ExpHeatMap2[y,x+1] <- length(which((ExpDisp2[,"x1"] == (ExpMin2 + x)) &
                                                  ExpDisp2[,"y1"] == y))
     }
}
image2D(t(ExpHeatMap2), x = ExpMin2:ExpMax2, y = 1:5, main = "Width = 10")


ExpMin3 <- min(ExpDisp3[,"x1"])
ExpMax3 <- max(ExpDisp3[,"x1"])
ExpHeatMap3 <- matrix(NA, nrow = 5, ncol = ExpMax3 - ExpMin3 + 1)
for(x in 0:(ExpMax3 - ExpMin3)){
     for(y in 1:5){
          ExpHeatMap3[y,x+1] <- length(which((ExpDisp3[,"x1"] == (ExpMin3 + x)) &
                                                  ExpDisp3[,"y1"] == y))
     }
}
image2D(t(ExpHeatMap3), x = ExpMin3:ExpMax3, y = 1:5, main = "Width = 20")

Check3 <- TRUE

if(Check1 & Check2 & Check3){
     writeLines("Pass", con = "ModelTesting/Disperse.txt")
} else{
     writeLines("Fail", con = "ModelTesting/Disperse.txt")
}

