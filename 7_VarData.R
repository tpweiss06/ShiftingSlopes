# This script will extract information on the per patch variability in fitness,
#    and the variability in mean trait values and abundance across the y dimension.

# CV = standard deviation / mean


# Set the number of processors available for this script
nProc <- 4*24

# Set the working directory
setwd("~/ShiftingSlopes/")
source("~/ShiftingSlopes/SimFunctions.R")
library(parallel)
library(Rmpi)
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")
NumSims <- 200

# Create a data frame to index each simulation block in the above arrays
TraitIndices <- expand.grid(params = 1:9, sim = 1:NumSims)

# Write a function to be passed to various nodes
VarExtract <- function(i){
     # Sort out the parameter combination and simulation under consideration
     Param <- TraitIndices$params[i]
     ParamDirectory <- paste("~/ShiftingSlopes/MainSim/Params", Param, "/", sep = "")
     AllSims <- list.files(ParamDirectory)
     SimID <- AllSims[TraitIndices$sim[i]]
     
     # Load the corresponding population matrix
     InFile <- paste("~/ShiftingSlopes/MainSim/Params", Param, "/", SimID,
                     "/SummaryStats.csv", sep = "")
     SimData <- read.csv(InFile)
     PopData <- subset(SimData, gen == 2000, select = c(x, y, abund, muFit, muDisp, sigmaFitPhen))
     return(PopData)
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, "TraitIndices")

# Run the simulations
SimVec <- 1:dim(TraitIndices)[1]
SimVar <- clusterApply(cl, x = SimVec, fun = VarExtract)

# For each scenario (1-9), I want a data frame with columns: x, abund, muFit, muDisp, sigmaFitPhen
# All of these columns except the last will be the variance across the y dimension for
#    these per patch quantities. sigmaFitPhen will be the mean across the y dimension
#    to get a sense of the average variance in reproductive success.

# Create the master list to hold all the values and the other necessary objects
VarList <- vector(length = 9, mode = "list")
for(p in 1:9){
     VarList[[p]] <- data.frame(x = rep(NA, 1), abund = rep(NA, 1), muFit = rep(NA, 1), muDisp = rep(NA, 1),
                                sigmaFitPhen = rep(NA, 1))
}

# Now populate the MasterList
ParamIndices <- rep(1, 9)
for(i in SimVec){
     CurParam <- TraitIndices$params[i]
     SimData <- SimVar[[i]]
     xVals <- unique(SimData$x)
     VarData <- data.frame(x = xVals, abund = rep(NA, length(xVals)), muFit = rep(NA, length(xVals)),
                           muDisp = rep(NA, length(xVals)), sigmaFitPhen = rep(NA, length(xVals)))
     for(j in 1:nrow(VarData)){
          xData <- subset(SimData, x == VarData$x[j])
          VarData$abund[j] <- sd(xData$abund, na.rm = TRUE) / mean(xData$abund, na.rm = TRUE)
          VarData$muFit[j] <- sd(xData$muFit, na.rm = TRUE) / mean(xData$muFit, na.rm = TRUE)
          VarData$muDisp[j] <- sd(xData$muDisp, na.rm = TRUE) / mean(xData$muDisp, na.rm = TRUE)
          VarData$sigmaFitPhen[j] <- mean(xData$sigmaFitPhen, na.rm = TRUE)
     }
     VarList[[CurParam]] <- rbind(VarList[[CurParam]], VarData)
}

# Save the output
save(VarList, file = "~/ShiftingSlopes/VarData.rdata")

