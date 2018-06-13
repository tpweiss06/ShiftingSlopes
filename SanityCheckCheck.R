# This script will double check that the sanity check simulations didn't go
#    extinct

SanitySims <- list.files("~/ShiftingSlopes/SanityCheck")
MaxGen <- rep(NA, length(SanitySims))
for(i in 1:length(SanitySims)){
     InFile <- paste("~/ShiftingSlopes/SanityCheck/", SanitySims[i], 
                     "/SummaryStats.csv", sep = "")
     SimData <- read.csv(InFile)
     MaxGen[i] <- max(SimData$gen)
}

MaxGen