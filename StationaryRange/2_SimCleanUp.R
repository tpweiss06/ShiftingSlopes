# Make a script to clean up unfinished simulations, extinct simulations, and too
#    many simulations

for(p in 1:9){
     ParamDirectory <- paste("~/ShiftingSlopes/StationaryRange/Params", p, "/", sep = "")
     SimFolders <- list.files(ParamDirectory)
     for(i in SimFolders){
          SimDirectory <- paste(ParamDirectory, i, sep = "")
          SimFiles <- list.files(SimDirectory)
          NumFiles <- length(SimFiles)
          if(NumFiles < 3){
               SysCommand1 <- paste("rm -r ", SimDirectory, sep = "")
               system(SysCommand1)
          } else{
               PopMatFile <- paste(SimDirectory, "/PopMat.csv")
               PopMat <- read.csv(PopMatFile)
               if(dim(PopMat)[1] == 0){
                    SysCommand2 <- paste("rm -r ", SimDirectory, sep = "")
                    system(SysCommand2)
               }
          }
     }
     NewSimFolders <- list.files(ParamDirectory)
     NumSims <- length(NewSimFolders)
     if(NumSims > 100){
          SimsToDelete <- NumSims - 100
          for(i in 1:SimsToDelete){
               SysCommand3 <- paste("rm -r ", ParamDirectory, NewSimFolders[i], sep = "")
               system(SysCommand3)
          }
     }
}