# This script will extract the realized fitness of individuals throughout their
#    ranges every 10 generations

# Set the number of processors available for this script 
nProc <- 2*24

# Set the working directory
setwd("~/ShiftingSlopes/")
source("~/ShiftingSlopes/SimFunctions.R")
library(parallel)
library(Rmpi)
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")
load("ExtSimIDs.rdata")

# Create the vectors to be passed to the function in the data extraction
SimVec <- 1:27
SpeedWords <- c("Slow", "MainSim", "Fast")
SpeedVec <- rep(SpeedWords, 9)
ParamVec <- rep(1:9, each = 3)

# Create other useful objects to pass to the nodes
GenSteps <- 10
LengthShift <- 100
GenSeq <- seq(0, LengthShift, by = GenSteps)

FitExtract <- function(i){
     # Set the current speed and parameter combination and isolate the relevant
     #    simulation IDs.
     SpeedWord <- SpeedVec[i]
     ParamCombo <- ParamVec[i]
     CurSims <- subset(SimIDs, Params == ParamCombo)
     
     # Source a representative parameter file for useful values
     ParamFile <- paste("~/ShiftingSlopes/", SpeedWord, "/Params", ParamCombo, 
                        "/", CurSims$SimID[1], "/parameters.R", sep = "")
     source(ParamFile)
     
     # Create the data frame for all simulations
     MeanFits <- data.frame(x = NA, g = NA, Mismatch = NA, lwr = NA, upr = NA, Zopt = NA)
     
     # Now loop through each row of the data frame to fill it in
     for(g in GenSeq){
          # Create a temporary data frame to hold the results from all simulations
          #    for each time point, which will then be condensed across simulations
          #    for the main data frames
          Temp <- data.frame(x = NA, Mismatch = NA, Zopt = NA)
          for(s in 1:nrow(CurSims)){
               if(g == 0){
                    InFile <- paste("~/ShiftingSlopes/MainSim/Params", ParamCombo,
                                    "/", CurSims$SimID[s], "/SummaryStats.csv", sep = "")
                    SumStats <- read.csv(InFile)
                    CurSumStats <- subset(SumStats, gen == 2000)
               } else{
                    InFile <- paste("~/ShiftingSlopes/", SpeedWord, "/Params", ParamCombo,
                                    "/", CurSims$SimID[s], "/SummaryStats.csv", sep = "")
                    SumStats <- read.csv(InFile)
                    CurSumStats <- subset(SumStats, gen == (g + BurnIn))
               }
               # Isolate the x coordinates from the current data and loop through
               #    them, adding their data to the temporary data frames
               xSeq <- unique(CurSumStats$x)
               for(j in xSeq){
                    LocalData <- subset(CurSumStats, x == j)
                    # Calculate the mismatch between optimum and genotype and
                    #    its direction (i.e. positive or negative)
                    Zopt <- lambda * (j * eta - LocalData$beta[1])
                    Maladaptation <- LocalData$muFit
                    # Add the simulation data to the Temp data frame
                    Temp <- rbind(c(j, mean(Maladaptation), Zopt), Temp)
               }
          }
          # remove the final row of NA's from Temp
          Temp <- Temp[-nrow(Temp),]
          # Find the unique x values, cycle through them, calculate mean and
          #    quantiles to go in master data frames
          Xvals <- unique(Temp$x)
          for(j in Xvals){
               TempData <- subset(Temp, x == j)
               FitQuants <- quantile(TempData$Mismatch, probs = c(0.25, 0.75))
               NewData <- c(j, g, mean(TempData$Mismatch), FitQuants, TempData$Zopt[1])
               MeanFits <- rbind(NewData, MeanFits)
          }
     }
     # Remove the final row of NA values
     Mismatch <- MeanFits[-nrow(MeanFits),]
     return(Mismatch)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("SpeedVec", "ParamVec", "SimIDs", "GenSeq"))

# Run the function on the cluster
SimFits <- clusterApply(cl, x = SimVec, fun = FitExtract)

# Process the results into an appropriate object for the later graphing script
MismatchData <- vector(mode = "list", length = 9)
for(p in 1:9){
     MismatchData[[p]] <- vector(mode = "list", length = 3)
     for(v in 1:3){
          CurIndex <- which((SpeedVec == SpeedWords[v]) & (ParamVec == p))
          CurResults <- SimFits[[CurIndex]]
          MismatchData[[p]][[v]] <- CurResults
     }
}

# Save the output
save(MismatchData, file = "~/ShiftingSlopes/MismatchData.rdata")

