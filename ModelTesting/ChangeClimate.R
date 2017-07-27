# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# To test the ChangeClimate function, first choose simple values of climate 
#    speed and initial conditions so that the answer is obvious from manual
#    calculations. Then supply a negative value of speed and ensure a helpful
#    error message is returned, followed by a length of time of 0 to ensure
#    that an empty vector is returned.

# First try a straitforward climate change
Beta1 <- 0
ClimateSpeed <- 1
ChangeDuration <- 10
# These values should produce a vector of length 10 with values 1, 2, 3, etc.

TestRun <- ChangeClimate(BetaInit = Beta1, LengthShift = ChangeDuration, 
                         ClimSpeed = ClimateSpeed)

if(sum(TestRun == 1:10) == 10){
     FirstCheck <- TRUE
} else{
     FirstCheck <- FALSE
}

# Next try a negative speed
TestNegative <- ChangeClimate(BetaInit = Beta1, LengthShift = ChangeDuration,
                              ClimSpeed = -1)
if(is.null(TestNegative)){
     SecondCheck <- TRUE
} else{
     SecondCheck <- FALSE
}

# Finally try a shift length of 0
Test0 <- ChangeClimate(BetaInit = Beta1, LengthShift = 0, ClimSpeed = ClimateSpeed)
if(is.null(Test0)){
     ThirdCheck <- TRUE
} else{
     ThirdCheck <- FALSE
}

if(FirstCheck & SecondCheck & ThirdCheck){
     writeLines("Pass", con = "ModelTesting/ChangeClimate.txt")
} else{
     writeLines("Fail", con = "ModelTesting/ChangeClimate.txt")
}