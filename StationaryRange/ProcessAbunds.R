# This script will process and save various abundance quantities (listed below) 
#    for subsequent graphing.
#    MeanAbunds: the mean abundance over all simulations in each individual patch
#                   at each time step.
#    SigmaAbunds: the standard deviation of abundances across simulations for each
#                   patch.
#    SectorMeans: the mean abundance over all simulations in each landscape sector 
#              (all patches in a column with the same range capacity value).
#    SectorSigma: the standard deviation of all simulations and patches within
#              a landscape sector.

# Set the working directory and declare which paramater combination to work with
setwd("~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/")
RangeParams <- 1

# Load in the simulation abundance data
load(paste("SimData/Params", RangeParams, "Abunds.rdata", sep = ""))
nWidth <- dim(Abunds[[1]])[1]
nLength <- dim(Abunds[[1]])[2]
nTime <- dim(Abunds[[1]])[3]

# Create the objects to hold the various quantities I need
MeanAbunds <- array(NA, dim = dim(Abunds[[1]]))
SigmaAbunds <- array(NA, dim = dim(Abunds[[1]]))
SectorMeans <- matrix(NA, nrow = nLength, ncol = nTime)
SectorSigma <- matrix(NA, nrow = nLength, ncol = nTime)

# Now make functions to extract the abundances of each patch in each simulation
GetAbunds <- function(wid, len, time, AbundArray){
     return(Abund)
}

# Now loop through all simulations and fill in these quantities


for(k in 1:nLength){
     for(l in 1:nTime){
          PooledAbunds <- NULL
          for(j in 1:nWidth){
               TempAbunds <- rep(NA, length(Abunds))
               for(i in 1:length(Abunds)){
                    TempAbunds[i] <- Abunds[[i]][j,k,l]
               }
               MeanAbunds[j,k,l] <- mean(TempAbunds)
               SigmaAbunds[j,k,l] <- sqrt(var(TempAbunds))
               PooledAbunds <- c(PooledAbunds, TempAbunds)
          }
          SectorMeans[k,l] <- mean(MeanAbunds[,k,l])
          SectorSigma[k,l] <- sqrt(var(PooledAbunds))
          rm(PooledAbunds)
     }
}

OutFile <- paste("SimData/Params", RangeParams, "ExtractedAbunds.rdata", sep = "")
save(MeanAbunds, SigmaAbunds, SectorMeans, SectorSigma, file = OutFile)