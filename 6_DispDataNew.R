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

DispExtract <- function(i){
     # Set the current speed and parameter combination and isolate the relevant
     #    simulation IDs.
     SpeedWord <- SpeedVec[i]
     ParamCombo <- ParamVec[i]
     CurSims <- subset(SimIDs, Params == ParamCombo)
     
     # Source a representative parameter file for useful values
     ParamFile <- paste("~/ShiftingSlopes/", SpeedWord, "/Params", ParamCombo, 
                        "/", CurSims$SimID[1], "/parameters.R", sep = "")
     source(ParamFile)
     
     # Create the data frames for extant and extinct simulations
     ExtantDisps <- data.frame(x = NA, g = NA, dBar = NA, lwr = NA, upr = NA)
     ExtinctDisps <- data.frame(x = NA, g = NA, dBar = NA, lwr = NA, upr = NA)
     
     # Now loop through each row of the data frame to fill it in
     for(g in GenSeq){
          # Create temporary data frames to hold the results from all simulations
          #    for each time point, which will then be condensed across simulations
          #    for the main data frames
          TempExtant <- data.frame(x = NA, dBar = NA)
          TempExtinct <- data.frame(x = NA, dBar = NA)
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
                    # Decide if this simulation goes in extinct or extant data
                    if(SpeedWord == "Slow"){
                         if(CurSims$Slow[s] == 1){
                              TempExtant <- rbind(c(j, mean(LocalData$muDisp)), TempExtant)
                         }else{
                              TempExtinct <- rbind(c(j, mean(LocalData$muDisp)), TempExtinct)
                         }
                    } else if(SpeedWord == "MainSim"){
                         if(CurSims$Moderate[s] == 1){
                              TempExtant <- rbind(c(j, mean(LocalData$muDisp)), TempExtant)
                         }else{
                              TempExtinct <- rbind(c(j, mean(LocalData$muDisp)), TempExtinct)
                         }
                    } else{
                         if(CurSims$Fast[s] == 1){
                              TempExtant <- rbind(c(j, mean(LocalData$muDisp)), TempExtant)
                         }else{
                              TempExtinct <- rbind(c(j, mean(LocalData$muDisp)), TempExtinct)
                         }
                    }
               }
          }
          # remove the final row of NA's from TempExtant and TempExtinct
          TempExtant <- TempExtant[-nrow(TempExtant),]
          TempExtinct <- TempExtinct[-nrow(TempExtinct),]
          # Find the unique x values, cycle through them, calculate mean and
          #    quantiles to go in master data frames
          ExtantXvals <- unique(TempExtant$x)
          ExtinctXvals <- unique(TempExtinct$x)
          if(nrow(TempExtant) > 0){
               for(j in ExtantXvals){
                    TempData <- subset(TempExtant, x == j)
                    DispQuants <- quantile(TempData$dBar, probs = c(0.25, 0.75))
                    NewData <- c(j, g, mean(TempData$dBar), DispQuants)
                    ExtantDisps <- rbind(NewData, ExtantDisps)
               }
          }
          if(nrow(TempExtinct) > 0){
               for(j in ExtinctXvals){
                    TempData <- subset(TempExtinct, x == j)
                    DispQuants <- quantile(TempData$dBar, probs = c(0.25, 0.75))
                    NewData <- c(j, g, mean(TempData$dBar), DispQuants)
                    ExtinctDisps <- rbind(NewData, ExtinctDisps)
               }
          }
     }
     # Remove the final row of NA values
     ExtantDisps <- ExtantDisps[-nrow(ExtantDisps),]
     ExtinctDisps <- ExtinctDisps[-nrow(ExtinctDisps),]
     # Combine and return the results
     Results <- list(Extant = ExtantDisps, Extinct = ExtinctDisps)
     return(Results)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("SpeedVec", "ParamVec", "SimIDs", "GenSeq"))

# Run the function on the cluster
SimDisps <- clusterApply(cl, x = SimVec, fun = DispExtract)

# Process the results into an appropriate object for the later graphing script
DispData <- vector(mode = "list", length = 9)
for(p in 1:9){
     DispData[[p]] <- vector(mode = "list", length = 3)
     for(v in 1:3){
          CurIndex <- which((SpeedVec == SpeedWords[v]) & (ParamVec == p))
          CurResults <- SimDisps[[CurIndex]]
          DispData[[p]][[v]] <- CurResults
     }
}

# Save the output
save(DispData, file = "~/ShiftingSlopes/DispDataNew.rdata")

