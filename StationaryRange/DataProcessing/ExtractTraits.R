# This script will extract the abundance information from the simulations

setwd("/home/shawa/cweissle/ShiftingSlopes/StationaryRange")
library(parallel)
library(Rmpi)
nProc <- 24

# Set the initial parameters for which range parameters we are using
RangeParams <- 1
SimFolder <- paste("Params", RangeParams, sep = "")
AllSimIDs <- list.files(SimFolder)
source(paste(SimFolder, AllSimIDs[1], "parameters.R", sep = "/"))

# Create objects to use in the extraction function below
TrialLength <- 1001
BetaPos <- TrialLength %/% 2

# Create a character vector with the names of all the objects needed on the nodes
ExportVec <- c("TrialLength", "BetaPos", "RangeParams", "width", "BurnIn")

# Create a function to extract the abundance information I want from each simulation
TraitExtract <- function(SimID){
     InFile <- paste("/home/shawa/cweissle/ShiftingSlopes/StationaryRange/Params",
                     RangeParams, "/", SimID, "/SummaryStats.csv", sep = "")
     SimData <- read.csv(InFile)
     MuFit <- array(0, dim = c(width, TrialLength, BurnIn))
     SigmaFit <- array(0, dim = c(width, TrialLength, BurnIn))
     MuDisp <- array(0, dim = c(width, TrialLength, BurnIn))
     SigmaDisp <- array(0, dim = c(width, TrialLength, BurnIn))
     for(g in 1:BurnIn){
          CurGen <- subset(SimData, (gen == g) & (abund > 0))
          xRange <- range(CurGen$x)
          xSeq <- seq(xRange[1], xRange[2], by = 1)
          for(j in xSeq){
               CurCol <- subset(CurGen, x == j)
               xArrInd <- BetaPos + j
               for(k in 1:width){
                    CurPatch <- subset(CurCol, y == k)
                    if(dim(CurPatch)[1] == 1){
                         MuFit[k, xArrInd, g] <- CurPatch$muFit
                         SigmaFit[k, xArrInd, g] <- CurPatch$sigmaFit
                         MuDisp[k, xArrInd, g] <- CurPatch$muDisp
                         SigmaDisp[k, xArrInd, g] <- CurPatch$sigmaDisp
                    }
               }
          }
     }
     Traits <- list(MuFit, SigmaFit, MuDisp, SigmaDisp)
     return(Traits)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, ExportVec)

# Run the simulations
SimTraits <- clusterApply(cl, x = AllSimIDs, fun = TraitExtract)

# Save the extracted abundance output
save(SimTraits, file = paste("Params", RangeParams, "Traits.rdata", sep = ""))

# The code will then shut down the cluster after this line.
