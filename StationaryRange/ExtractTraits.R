# This script will extract the trait information from the simulations

setwd("/home/shawa/cweissle/ShiftingSlopes/StationaryRange")

# Set the initial parameters for which range parameters we are using
RangeParams <- 1
SimFolder <- paste("Params", RangeParams, sep = "")
AllSimIDs <- list.files(SimFolder)
source(paste(SimFolder, AllSimIDs[1], "parameters.R", sep = "/"))

# Once the abunds script works without errors, use that as a template here to also
#    extract the trait data.