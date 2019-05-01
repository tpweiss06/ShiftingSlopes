# This script will extract the dispersal phenotypes of individuals throughout their
#    ranges every 10 generations

# Set the number of processors available for this script 
nProc <- 2*24

# Set the working directory
setwd("~/ShiftingSlopes/")
source("~/ShiftingSlopes/SimFunctions.R")
library(parallel)
library(Rmpi)
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")

# Create the vectors to be passed to the function in the data extraction
SimVec <- 1:9

# Create other useful objects to pass to the nodes
GenSteps <- 10
LengthShift <- 2000
GenSeq <- seq(10, LengthShift, by = GenSteps)

DispExtract <- function(i){
     # Set the current parameter combination and isolate the relevant
     #    simulation IDs.
     Prefix <- paste("~/ShiftingSlopes/EquilibriumTest/Params", i, sep = "")
     CurSims <- list.files(Prefix)
     
     # Source a representative parameter file for useful values
     ParamFile <- paste(Prefix, CurSims[1], "parameters.R", sep = "/")
     source(ParamFile)
     
     # Create the data frame for all simulations
     MeanDisps <- data.frame(x = NA, g = NA, dBar = NA, lwr = NA, upr = NA)
     
     # Now loop through each row of the data frame to fill it in
     for(g in GenSeq){
          # Create a temporary data frame to hold the results from all simulations
          #    for each time point, which will then be condensed across simulations
          #    for the main data frames
          Temp <- data.frame(x = NA, dBar = NA)
          for(s in 1:length(CurSims)){
               InFile <- paste(Prefix, CurSims[s], "SummaryStats.csv", sep = "/")
               SumStats <- read.csv(InFile)
               CurSumStats <- subset(SumStats, gen == g)
               # Isolate the x coordinates from the current data and loop through
               #    them, adding their data to the temporary data frames
               xSeq <- unique(CurSumStats$x)
               for(j in xSeq){
                    LocalData <- subset(CurSumStats, x == j)
                    # Add the simulation data to the Temp data frame
                    Temp <- rbind(c(j, mean(LocalData$muDisp)), Temp)
               }
          }
          # remove the final row of NA's from Temp
          Temp <- Temp[-nrow(Temp),]
          # Find the unique x values, cycle through them, calculate mean and
          #    quantiles to go in master data frames
          Xvals <- unique(Temp$x)
          for(j in Xvals){
               TempData <- subset(Temp, x == j)
               DispQuants <- quantile(TempData$dBar, probs = c(0.25, 0.75))
               NewData <- c(j, g, mean(TempData$dBar), DispQuants)
               MeanDisps <- rbind(NewData, MeanDisps)
          }
     }
     # Remove the final row of NA values
     MeanDisps <- MeanDisps[-nrow(MeanDisps),]
     return(MeanDisps)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, "GenSeq")

# Run the function on the cluster
DispData <- clusterApply(cl, x = SimVec, fun = DispExtract)

# Save the output
save(DispData, file = "~/ShiftingSlopes/EquilibriumDispData.rdata")

