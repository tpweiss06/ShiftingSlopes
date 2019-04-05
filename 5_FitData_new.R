# This script will extract the realized fitness of individuals throughout their
#    ranges in the generation immediately prior to climate change

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
xSeq <- -40:40

FitExtract <- function(i){
     # Set the current speed and parameter combination and isolate the relevant
     #    simulation IDs
     SpeedWord <- SpeedVec[i]
     ParamCombo <- ParamVec[i]
     CurSims <- subset(SimIDs, Params == ParamCombo)
     
     # Source a representative parameter file for useful values
     ParamFile <- paste("~/ShiftingSlopes/MainSim/Params", ParamCombo, "/", CurSims$SimID[1],
                        "/parameters.R", sep = "")
     source(ParamFile)
     
     # Create the data frames for extant and extinct simulations
     ExtantFits <- data.frame(x = xSeq, wBar = rep(NA, length(xSeq)),
                              lwr = rep(NA, length(xSeq)), upr = rep(NA, length(xSeq)))
     ExtinctFits <- data.frame(x = xSeq, wBar = rep(NA, length(xSeq)),
                              lwr = rep(NA, length(xSeq)), upr = rep(NA, length(xSeq)))
     
     # Now loop through each row of the data frame to fill it in
     for(j in 1:length(xSeq)){
          ExtantVals <- NULL
          ExtinctVals <- NULL
          for(s in 1:nrow(CurSims)){
               FilePath <- paste("~/ShiftingSlopes/", SpeedWord, "/Params", ParamCombo, "/", CurSims$SimID[s], "/InitialPopMat.csv", sep = "")
               SimData <- read.csv(InFile)
               if( (j == 1) & (s == 1)){ # Only need to do this once
                    FitCols <- grep(pattern = "fit", x = names(SimData)) 
               }
               CurXvals <- subset(SimData, x0 == xSeq[j])
               if(dim(CurXvals) > 0){
                    Zopt <- lambda * (CurXvals$x0 * eta - BetaInit)
                    EnvNiche <- rowSums(CurXvals[, FitCols])
                    RelFits <- exp(-1 * (EnvNiche - Zopt)^2 / (2*omega^2))
                    if(SpeedWord == "Slow"){
                         if(CurSims$Slow[s] == 1){
                              ExtantVals <- c(ExtantVals, RelFits)
                         }else{
                              ExtinctVals <- c(ExtinctVals, RelFits)
                         }
                    } else if(SpeedWord = "MainSim"){
                         if(CurSims$Moderate[s] == 1){
                              ExtantVals <- c(ExtantVals, RelFits)
                         }else{
                              ExtinctVals <- c(ExtinctVals, RelFits)
                         }
                    } else{
                         if(CurSims$Fast[s] == 1){
                              ExtantVals <- c(ExtantVals, RelFits)
                         }else{
                              ExtinctVals <- c(ExtinctVals, RelFits)
                         }
                    }
               }
          }
          if(!is.null(ExtantVals)){
               ExtantFits$wBar[j] <- mean(ExtantVals)
               Quantiles <- quantile(x = ExtantVals, probs = c(0.25, 0.75))
               ExtantFits$lwr[j] <- Quantiles[1]
               ExtantFits$upr[j] <- Quantiles[2]
          }
          if(!is.null(ExtinctVals)){
               ExtinctFits$wBar[j] <- mean(ExtinctVals)
               Quantiles <- quantile(x = ExtinctVals, probs = c(0.25, 0.75))
               ExtinctFits$lwr[j] <- Quantiles[1]
               ExtinctFits$upr[j] <- Quantiles[2]
          }
     }
     Results <- list(Extant = ExtantFits, Extinct = ExtinctFits)
     return(Results)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("SpeedVec", "ParamVec", "SimIDs", "xSeq"))

# Run the function on the cluster
SimFits <- clusterApply(cl, x = SimVec, fun = FitExtract)

# Process the results into an appropriate object for the later graphing script
FitData <- vector(mode = "list", length = 9)
for(p in 1:9){
     FitData[[p]] <- vector(mode = "list", length = 3)
     for(v in 1:3){
          CurIndex <- which((SpeedVec == SpeedWords[v]) & (ParamVec == p))
          CurResults <- SimFits[[CurIndex]]
          FitData[[p]][[v]] <- CurResults
     }
}

# Save the output
save(FitData, file = "~/ShiftingSlopes/FitDataNew.rdata")

