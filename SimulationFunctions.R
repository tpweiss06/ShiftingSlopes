# This script contains all of the functions necessary to run a simulation of
#    population evolution and movement in response to climate change. The 
#    functions are organized into three different categories: biological,
#    bookkeeping, and environmental. Each function is documented with a brief
#    description of the purpose of the function and a list and definition of 
#    the necessary function arguments before the function definition itself.
#    Each broad category of functions is delineated by the dashed section
#    headers and each individual function begins with a single line of 6 pound
#    signs (######)


# ------------------------------------------------------------------------------
# ---------------------------- Biological Functions ----------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# -------------------------- Environmental Functions ---------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# --------------------------- Bookkeeping Functions ----------------------------
# ------------------------------------------------------------------------------

######
# This function generates the matrices used to track individuals in the model
### INPUTS
# PopSize: a single number indicating the number of individuals the matrix needs
#    to track
# nFit: the number of loci defining an individual's fitness
# nDisp: the number of loci defining an individual's dispersal ability
# example: a boolean value indicating whether the function should simply return
#    an example character vector containing the column names and order used
### OUTPUTS
# This function creates a single, empty matrix wth PopSize rows and nFit + nDisp
#    + 5 columns with appropriate names (see function definition for details). 
#    If example = TRUE, however, the function simply returns a simple character
#    vector with the column names and layout.
MakePopMat <- function(PopSize, nFit, nDisp, example = FALSE){
     if(example){
          MatNames <- c("x0", "y0", "x1", "y1", "sex", "fit1", "...", "fitN",
                        "disp1", "...", "dispN")
          return(MatNames)
     }
     # First determine the number of columns and make the population matrix
     MatCols <- (nFit + nDisp + 5)
     PopMat <- matrix(NA, nrow = PopSize, ncol = MatCols)
     # Next create a vector of column names
     MatNames <- rep(NA, MatCols)
     MatNames[1:5] <- c("x0", "y0", "x1", "y1", "sex")
     for(n in 1:nFit){
          ColName <- paste("fit", n, sep = "")
          MatNames[5 + n] <- ColName
     }
     for(n in 1:nDisp){
          ColName <- paste("disp", n, sep = "")
          MatNames[5 + nFit + n] <- ColName
     }
     
     # Apply the names and return the matrix
     colnames(PopMat) <- MatNames
     return(PopMat)
}

######
# This function looks through the current working directory and finds a safe
#    name to use that will not overwrite anything else when saving simulation
#    results
### INPUTS
# ParentDirectory: the file path to the current working directory
# parallel: A boolean variable indicating whether the current simulations are
#    taking place over multiple nodes so that file names should include the 
#    name of the current node
### OUTPUTS
# The full file path for a new directory that the current simulation results will
#    be saved to
get_safe_ID <- function(ParentDirectory, parallel = FALSE){
     # Set the working directory to easily scan for directory names
     setwd(ParentDirectory)
     
     # Start by trying out 1 and then increase as necessary
     newID <- 1
     
     # If the simulations are being done in parallel, include the name of the
     #    current node in the new directory name
     if(parallel){
          NodeName <- Sys.getpid()
          DirName <- paste(NodeName, "_", "Sim", newID, sep = "")
          # Until we find a directory name that does not already exist, continue
          #    to increase the ID variable
          while(dir.exists(DirName)){
               newID <- newID + 1
               DirName <- paste(NodeName, "_", "Sim", newID, sep = "")
          }
     } else{
          DirName <- paste('Sim', newID, sep='')
          while(dir.exists(DirName)){
               newID <- newID + 1
               DirName <- paste('Sim', newID, sep='')
          }
     }
     
     # Now paste together the ParentDirectory with the new directory for the
     #    full file path
     FullPath <- paste(ParentDirectory, DirName, sep = "/")
     return(FullPath)
}

