# This script will extract the realized fitness of individuals throughout their
#    ranges in the generation immediately prior to climate change

# Set the number of processors available for this scrip
nProc <- 4*24

# Set the working directory
setwd("~/ShiftingSlopes/")
source("~/ShiftingSlopes/SimFunctions.R")
library(parallel)
library(Rmpi)
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")

# Create the master list to hold all the values and the other necessary objects
NumSims <- 200
MaxGen <- c(151, 2151, 151)
FilePaths <- c("~/ShiftingSlopes/Slow/Params", "~/ShiftingSlopes/MainSim/Params",
               "~/ShiftingSlopes/Fast/Params")
FitList <- vector(length = 9, mode = "list")
for(p in 1:9){
     FitList[[p]] <- vector(length = NumSims, mode = "list")
}

# Create a data frame to index each simulation block in the above arrays
TraitIndices <- expand.grid(params = 1:9, sim = 1:NumSims)

# Write a function to be passed to various nodes
FitExtract <- function(i){
     # Sort out the parameter combination and simulation under consideration
     Param <- TraitIndices$params[i]
     ParamDirectory <- paste("~/ShiftingSlopes/MainSim/Params", Param, "/", sep = "")
     AllSims <- list.files(ParamDirectory)
     SimID <- AllSims[TraitIndices$sim[i]]
     
     # Load the corresponding population matrix
     InFile <- paste("~/ShiftingSlopes/MainSim/Params", Param, "/", SimID,
                     "/InitialPopMat.csv", sep = "")
     ParamFile <- paste("~/ShiftingSlopes/MainSim/Params", Param, "/", SimID,
                        "/parameters.R", sep = "")
     SimData <- read.csv(InFile)
     source(ParamFile)
     
     # Convert the phenotype values to relative fitness values based on where they are
     #    in the range
     FitCols <- grep(pattern = "fit", x = names(SimData))
     Zopt <- lambda * (SimData$x0 * eta - BetaInit)
     EnvNiche <- rowSums(SimData[, FitCols])
     RelFits <- exp(-1 * (EnvNiche - Zopt)^2 / (2*omega^2))
     
     # Put it all together
     w <- matrix(NA, nrow = length(RelFits), ncol = 2)
     w[,1] <- SimData$x0
     w[,2] <- RelFits
     
     # Check if the simulation went extinct at each speed of climate change
     Ext <- rep(NA, 3)
     for(v in 1:3){
          SumFile <- paste(FilePaths[v], Param, "/", SimID, "/SummaryStats.csv", sep = "")
          SumStats <- read.csv(SumFile)
          GenExtinct <- max(SumStats$gen) + 1
          Ext[v] <- ifelse(GenExtinct < (MaxGen[v]), 0, 1)
     }
     return(list(w = w, Ext = Ext))
}

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("TraitIndices", "MaxGen", "FilePaths"))

# Run the simulations
SimVec <- 1:dim(TraitIndices)[1]
SimFit <- clusterApply(cl, x = SimVec, fun = FitExtract)

# Now populate the MasterList
ParamIndices <- rep(1, 9)
for(i in SimVec){
     CurParam <- TraitIndices$params[i]
     CurParamIndex <- ParamIndices[CurParam]
     FitList[[CurParam]][[CurParamIndex]] <- SimFit[[i]]
     ParamIndices[CurParam] <- ParamIndices[CurParam] + 1
}

# Save the output
save(FitList, file = "~/ShiftingSlopes/FitData.rdata")

