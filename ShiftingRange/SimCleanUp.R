# This script will walk through all the parameter folders and remove any
#    simulation subfolders that are incomplete (i.e. don't have all three
#    output files)

setwd("~/ShiftingSlopes/ShiftingRange/Slow")

for(p in 1:9){
     NewDirectory <- paste("Params", p, sep = "")
     SimFiles <- list.files(NewDirectory)
     for(j in 1:100){
          SubDirectory <- paste(NewDirectory, SimFiles[j], sep = "/")
          if(NumFiles < 3){
               SysCommand <- paste("rm -r ", SubDirectory, sep = "")
               system(SysCommand)
          }
     }
}

setwd("~/ShiftingSlopes/ShiftingRange/Med")

for(p in 1:9){
     NewDirectory <- paste("Params", p, sep = "")
     SimFiles <- list.files(NewDirectory)
     for(j in 1:100){
          SubDirectory <- paste(NewDirectory, SimFiles[j], sep = "/")
          if(NumFiles < 3){
               SysCommand <- paste("rm -r ", SubDirectory, sep = "")
               system(SysCommand)
          }
     }
}

setwd("~/ShiftingSlopes/ShiftingRange/Fast")

for(p in 1:9){
     NewDirectory <- paste("Params", p, sep = "")
     SimFiles <- list.files(NewDirectory)
     for(j in 1:100){
          SubDirectory <- paste(NewDirectory, SimFiles[j], sep = "/")
          if(NumFiles < 3){
               SysCommand <- paste("rm -r ", SubDirectory, sep = "")
               system(SysCommand)
          }
     }
}
