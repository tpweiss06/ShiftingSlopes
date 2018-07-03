# This script will check for extinction before the period of climate change,
#    delete those simulations, and output the current number of simulations
#    that make it to equilibriation.

setwd("~/ShiftingSlopes/MainSim/")
NumSims <- rep(NA, 9)

for(p in 1:9){
     ParamDirectory <- paste("~/ShiftingSlopes/MainSim/Params", p, sep = "")
     ParamSims <- list.files(ParamDirectory)
     for(i in 1:length(ParamSims)){
          NumFiles <- length(list.files(paste(ParamDirectory, ParamSims[i], sep = "/")))
          if(NumFiles != 4){
               SysCommand <- paste("rm -r ", ParamDirectory, "/", ParamSims[i], sep = "")
               system(SysCommand)
          } else{
               InFile <- paste(ParamDirectory, ParamSims[i], "InitialPopMat.csv", sep = "/")
               PopMat <- read.csv(InFile)
               if(dim(PopMat)[1] == 0){
                    SysCommand <- paste("rm -r ", ParamDirectory, "/", ParamSims[i], sep = "")
                    system(SysCommand)
               }
          }
     }
     NumSims[p] <- length(list.files(ParamDirectory))
}

NumSims