# This script will walk through all the parameter folders and remove any
#    simulation subfolders that are incomplete (i.e. don't have all three
#    output files)

SpeedWords <- c("Slow", "Med", "Fast")
NumSims <- matrix(NA, nrow = 3, ncol = 9)
for(s in 1:3){
     for(p in 1:9){
          ParamDirectory <- paste("~/ShiftingSlopes/ShiftingRange/", SpeedWords[s],
                                  "/Params", p, "/", sep = "")
          SimFiles <- list.files(ParamDirectory)
          for(i in SimFiles){
               SimDirectory <- paste(ParamDirectory, i, sep = "")
               NumFiles <- length(list.files(SimDirectory))
               if(NumFiles < 3){
                    SysCommand <- paste("rm -r ", SimDirectory, sep = "")
                    system(SysCommand)
               }
          }
          RemainingSims <- list.files(ParamDirectory)
          NumSims[s,p] <- length(RemainingSims)
     }
}

NumSims

100 - NumSims

