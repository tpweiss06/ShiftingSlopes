# This script will walk through all the parameter folders and remove any
#    simulation subfolders that are incomplete (i.e. don't have all three
#    output files)

setwd("~/ShiftingSlopes/ShiftingRange/Slow")

for(p in 1:4){
     NewDirectory <- paste("Params", p, sep = "")
     #print(NewDirectory)
     SimFiles <- list.files(NewDirectory)
     for(j in 1:100){
          SubDirectory <- paste(NewDirectory, SimFiles[j], sep = "/")
          #NumFiles <- length(list.files(SubDirectory))
          #if(NumFiles < 3){
               SysCommand <- paste("rm -r ", SubDirectory, sep = "")
               system(SysCommand)
     	  #     print(SysCommand)	
          #}
     }
}

#NewDirectory <- "Params5"
#SimFiles <- list.files(NewDirectory)
#for(j in 1:62){
#     SubDirectory <- paste(NewDirectory, SimFiles[j], sep = "/")
#     SysCommand <- paste("rm -r ", SubDirectory, sep = "")
#     system(SysCommand)
#}
