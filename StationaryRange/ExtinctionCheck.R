# Make a script to test if a simulation went extinct in the stationary scenario
#    and then remove that simulation

setwd("~/ShiftingSlopes/StationaryRange")

NumExtinct <- rep(NA, 9)
for(p in 1:9){
     NewDirectory <- paste("Params", p, sep = "")
     SimFiles <- list.files(NewDirectory)
     Extinct <- rep(NA, 100)
     for(j in 1:100){
          PopMatFile <- InFile <- paste("~/ShiftingSlopes/StationaryRange/Params", p, "/", 
                                    SimFiles[j], "/PopMat.csv", sep = "")
          PopMat <- read.csv(PopMatFile)
          if(dim(PopMat)[1] == 0){
               Extinct[j] <- TRUE
               SysCommand <- paste("rm -r ~/ShiftingSlopes/StationaryRange/Params", p, "/",
                                   SimFiles[j], sep = "")
               system(SysCommand)
          } else{
               Extinct[j] <- FALSE
          }
     }
     NumExtinct[p] <- sum(Extinct)
}

print(NumExtinct)
