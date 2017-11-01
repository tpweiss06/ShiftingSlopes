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
AbundExtract <- function(SimID){
     InFile <- paste("/home/shawa/cweissle/ShiftingSlopes/StationaryRange/Params",
                     RangeParams, "/", SimID, "/SummaryStats.csv", sep = "")
     SimData <- read.csv(InFile)
     Abunds <- array(0, dim = c(width, TrialLength, BurnIn))
     for(g in 1:BurnIn){
          CurGen <- subset(SimData, gen == g)
          for(y in 1:width){
               CurRow <- subset(CurGen, y == y)
               CurCols <- CurRow$x
               for(x in 1:length(CurCols)){
                    PopSize <- CurRow$abund[x]
                    xArrInd <- BetaPos + CurCols[x]
                    if((xArrInd > 0) & (xArrInd <= TrialLength)){
                         Abunds[y, xArrInd, g]
                    }
               }
          }
     }
     return(Abunds)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, ExportVec)

# Run the simulations
SimAbunds <- clusterApply(cl, x = AllSimIDs, fun = AbundExtract)

# Save the extracted abundance output
save(SimAbunds, paste("Params", RangeParams, "Abunds.rdata", sep = ""))

# The code will then shut down the cluster after this line.