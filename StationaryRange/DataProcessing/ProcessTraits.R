# This script will process and save various trait quantities (listed below) 
#    for subsequent graphing.
#    MeanMuFit: the mean MuFit value over all simulations.
#    MeanSigFit: the mean SigmaFit value over all simulations.
#    MeanMuDisp: the same quantity as above but for dispersal.
#    MeanSigDisp: the same quantity as above but for dispersal.
#    SectorMuFit: the mean MuFit over all patches in a landscape sector.
#    SectorSigFit: the mean SigmaFit over all patches in a landscape sector.
#    SectorMuDisp: the same quantity as above but for dispersal.
#    SectorSigDisp: the same quantity as above but for dispersal.

# Set the working directory and declare which paramater combination to work with
setwd("/home/shawa/cweissle/ShiftingSlopes/StationaryRange")
#setwd("~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/")
RangeParams <- 1

# Load in the simulation abundance data
#load(paste("SimData/Params", RangeParams, "/Params", RangeParams, "Traits.rdata", sep = ""))
load(paste("Params", RangeParams, "Traits.rdata", sep = ""))
nWidth <- dim(SimTraits[[1]][[1]])[1]
nLength <- dim(SimTraits[[1]][[1]])[2]
nTime <- dim(SimTraits[[1]][[1]])[3]

# Create the objects to hold the various quantities I need
MeanMuFit <- array(NA, dim = dim(SimTraits[[1]][[1]]))
MeanSigFit <- array(NA, dim = dim(SimTraits[[1]][[1]]))
MeanMuDisp <- array(NA, dim = dim(SimTraits[[1]][[1]]))
MeanSigDisp <- array(NA, dim = dim(SimTraits[[1]][[1]]))
SectorMuFit <- matrix(NA, nrow = nLength, ncol = nTime)
SectorSigFit <- matrix(NA, nrow = nLength, ncol = nTime)
SectorMuDisp <- matrix(NA, nrow = nLength, ncol = nTime)
SectorSigDisp <- matrix(NA, nrow = nLength, ncol = nTime)

# Now loop through all simulations and fill in these quantities
for(k in 1:nLength){
     for(l in 1:nTime){
          for(j in 1:nWidth){
               TempMuFit <- rep(NA, length(SimTraits))
               TempSigFit <- rep(NA, length(SimTraits))
               TempMuDisp <- rep(NA, length(SimTraits))
               TempSigDisp <- rep(NA, length(SimTraits))
               for(i in 1:length(SimTraits)){
                    TempMuFit[i] <- SimTraits[[i]][[1]][j,k,l]
                    TempSigFit[i] <- SimTraits[[i]][[2]][j,k,l]
                    TempMuDisp[i] <- SimTraits[[i]][[3]][j,k,l]
                    TempSigDisp[i] <- SimTraits[[i]][[4]][j,k,l]
               }
               MeanMuFit[j,k,l] <- mean(TempMuFit, na.rm = TRUE)
               MeanSigFit[j,k,l] <- mean(TempSigFit, na.rm = TRUE)
               MeanMuDisp[j,k,l] <- mean(TempMuDisp, na.rm = TRUE)
               MeanSigDisp[j,k,l] <- mean(TempSigDisp, na.rm = TRUE)
          }
          SectorMuFit[k,l] <- mean(MeanMuFit[,k,l], na.rm = TRUE)
          SectorSigFit[k,l] <- mean(MeanSigFit[,k,l], na.rm = TRUE)
          SectorMuDisp[k,l] <- mean(MeanMuDisp[,k,l], na.rm = TRUE)
          SectorSigDisp[k,l] <- mean(MeanSigDisp[,k,l], na.rm = TRUE)
     }
     #print(k / nLength)
}

#OutFile <- paste("SimData/Params", RangeParams, "/ExtractedTraits.rdata", sep = "")
OutFile <- paste("Params", RangeParams, "ExtractedTraits.rdata", sep = "")
save(MeanMuFit, MeanSigFit, MeanMuDisp, MeanSigDisp, SectorMuFit, SectorSigFit,
     SectorMuDisp, SectorSigDisp, file = OutFile)