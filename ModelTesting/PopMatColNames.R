# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# To test the PopMatColNames function, first ensure that the example functionality
#    works regardless of the other arguments. Then, check that it returns the
#    expected output for different combinations of number of alleles and
#    life history syndromes (monoecious vs. dioecious).

Test1 <- PopMatColNames(example = TRUE)
Test2 <- PopMatColNames(nFit = 10, nDisp, monoecious = FALSE, example = TRUE)
if(sum(Test1 == Test2) == length(Test1)){
     Check1 <- TRUE
} else{
     Check1 <- FALSE
}

Test3 <- PopMatColNames(nFit = 3, nDisp = 4, monoecious = TRUE)       # Should be 18
Test4 <- PopMatColNames(nFit = 1, nDisp = 1, monoecious = FALSE)      # Should be 9
Test5 <- PopMatColNames(nFit = 20, nDisp = 100, monoecious = TRUE)    # Should be 244

if( (length(Test3) == 18) & (length(Test4) == 9) & (length(Test5) == 244) ){
     Check2 <- TRUE
} else{
     Check2 <- FALSE
}     
     
if( Check1 & Check2 ){     
     writeLines("Pass", con = "ModelTesting/PopMatColNames.txt")
} else{
     writeLines("Fail", con = "ModelTesting/PopMatColNames.txt")
}
