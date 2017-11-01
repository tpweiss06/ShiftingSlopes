# This script will create some preliminary graphs for the stationary range
#    explorations.
# Topher Weiss-Lehman
# October 23, 2017

setwd("~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/")
library(plot3D)

# Set the initial parameters for which range parameters we are using
RangeParams <- 1
SimFolder <- paste("Params", RangeParams, sep = "")
AllSimIDs <- list.files(SimFolder)
source(paste(SimFolder, AllSimIDs[1], "parameters.R", sep = "/"))

# Create objects to hold all the abundance data from all simulations
TrialLength <- 1001
BetaPos <- TrialLength %/% 2
AllAbunds <- array(0, dim = c(width, TrialLength, BurnIn, length(AllSimIDs)))

# Run through all the simulations and store the data
for(i in 1:AllSimIDs){
     SimData <- read.csv(paste(SimFolder, "SummaryStats.csv", sep = "/"))
     for(g in 1:BurnIn){
          CurGen <- subset(SimData, gen == g)
          for(y in 1:width){
               CurRow <- subset(CurGen, y == y)
               CurCols <- CurRow$x
               for(x in 1:length(CurCols)){
                    PopSize <- CurRow$abund[x]
                    xArrInd <- BetaPos + CurCols[x]
                    if((xArrInd > 0) & (xArrInd <= TrialLength)){
                         AllAbunds[y, xArrInd, g, i]
                    }
               }
          }
     }
}

# Now create and populate objects to hold summary statistics
MeanAbunds <- array(0, dim = c(width, TrialLength, BurnIn))
StDevAbunds <-  array(0, dim = c(width, TrialLength, BurnIn))

for(w in 1:width){
     for(l in 1:TrialLength){
          for(g in 1:BurnIn){
               MeanAbunds[w,l,g] <- mean(AllAbunds[w, l, g,])
               StDevAbunds[w,l,g] <- sd(AllAbunds[w, l, g,])
          }
     }
}


# What I need: two arrays with dimensions of width, length, and time; one for 
#    the mean abundance and one for the standard deviation in abundance.
