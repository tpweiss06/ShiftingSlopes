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

###### CalcDispTrait
# This function will calculate a dispersal trait value for a set of individuals
#    based on their dispersal alleles.
### INPUTS
# population:  A vector of the population (i.e. row numbers) for which to 
#                   calculate a dispersal trait.
# PopMat:      The population matrix to use
# nDisp:       The number of dispersal loci in the current simulation
### OUTPUTS
# The function will create a matrix with two columns. The first column will 
#    correspond to the ID (row number) of the individual and the second column
#    will be the dispersal trait value. If the dispersal trait is calculated for
#    the entire population, as will normally be the case, then the first column
#    will simply correspond to the row number, but including this column allows
#    for the calculation of dispersal for a subset of the current population.
CalcDispTrait <- function(population, PopMat, nDisp){
     # First create an empty matrix for the dispersal trait values and put
     #    the population IDs in the first column
     DispTrait <- matrix(NA, nrow = length(population), ncol = 2)
     DispTrait[,1] <- population
     # Calculate which columns of PopMat correspond to dispersal. Since it is
     #    the last trait in the matrix, we don't need to know how many fitness
     #    traits occured before it, we simply subtract from the last column.
     DispColumns <- ncol(PopMat) - nDisp + 1
     # Now calculate each individual's dispersal loci sum
     LociSum <- rowSums(PopMat[population, DispColumns])
     # Finally, exponentiate this sum to ensure that the dispersal trait 
     #    (diffusion coefficient) can't become negative. There is no need to
     #    impose an upper limit because the range limit itself will impose a
     #    cost on dispersing too far.
     DispTrait[,2] <- exp(LociSum)
     return(DispTrait)
}

###### Disperse
# This function will use the previously calculated diffusion coefficients in a
#    Poisson process to perform the actual dispersal events.
### INPUTS
# PopMat:      The population matrix to use
# DispTrait:   The two column matrix created by the CalcDispTrait function
# width:       The number of discrete cells defining the width of the word
#                   being simulated
# DispTime:    The total time allowed for dispersal during the Poisson Process.
#                   This will be set to 1 by default because it really only
#                   adjusts the scale of the diffusion coefficients so there's
#                   no good reason to set to anything else unless we are trying
#                   to match scales with a given system or a different model.
### OUTPUTS
# The function will perform dispersal for all the individuals involved via a 
#    Poisson process in which each event is defined by an individual randomly
#    shifting in space one cell in any direction (i.e. unbiased). This function
#    assumes wrapping borders for the width direction to avoid edge effects.
Disperse <- function(PopMat, DispTrait, width, DispTime = 1){
     # First create the initial vector of event times based on each individual's
     #    diffusion coefficient.
     EventTimes <- rexp(n = nrow(DispTrait), rate = 1) / DispTrait[,2]
     # Find the first event and record the individual
     CurTime <- min(EventTimes)
     CurIndividual <- which(EventTimes == CurTime)
     # Record the current x and y coordinates for each individual prior to 
     #    dispersal
     CurX <- PopMat[DispTrait[,1], "x0"]
     CurY <- PopMat[DispTrait[,1], "y0"]
     # Now disperse individuals according to their event times until we exceed
     #    the allowed dispersal time
     StepDists <- c(-1, 1)
     while(CurTime <= DispTime){
          # Randomly choose a direction for dispersal
          XorY <- rbinom(n = 1, size = 1, prob = 0.5)
          if(XorY == 0){
               # Dispersal occurs in the x direction
               CurX[CurIndividual] <- CurX[CurIndividual] + sample(StepDists, size = 1)
          }
          if(XorY == 1){
               # Dispersal occurs in the y direction
               CurY[CurIndividual] <- CurY[CurIndividual] + sample(StepDists, size = 1)
          }
          # Add in the next event time for the current individual and reset
          #    the current time and individual
          EventTimes[CurIndividual] <- rexp(n = 1, rate = 1) / DispTrait[CurIndividual,2]
          CurTime <- min(EventTimes)
          CurIndividual <- which(EventTimes == CurTime)
     }
     # Use the width of the landscape to correct any x coordinates outside
     #    of the allowed boundaries
     while(min(CurX) < 0){
          CurX <- ifelse(CurX < 0, CurX + width, CurX)
     }
     CurX <- ifelse(CurX > width, CurX %% width, CurX)
     CurX <- ifelse(CurX == 0, width, CurX)
     
     # Finally, update the appropriate entries of the PopMat matrix and return
     #    the new version.
     PopMat[DispTrait[,1], "x1"] <- CurX
     PopMat[DispTrait[,1], "y1"] <- CurY
     return(PopMat)
}





# ------------------------------------------------------------------------------
# -------------------------- Environmental Functions ---------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# --------------------------- Bookkeeping Functions ----------------------------
# ------------------------------------------------------------------------------

###### MakePopMat
# This function generates the matrices used to track individuals in the model
### INPUTS
# PopSize:     a single number indicating the number of individuals the matrix 
#                   needs to track
# nFit:        the number of loci defining an individual's fitness
# nDisp:       the number of loci defining an individual's dispersal ability
# example:     a boolean value indicating whether the function should simply 
#                   return an example character vector containing the column 
#                   names and order used
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

###### get_safe_ID
# This function looks through the current working directory and finds a safe
#    name to use that will not overwrite anything else when saving simulation
#    results
### INPUTS
# ParentDirectory:  the file path to the current working directory
# parallel:         A boolean variable indicating whether the current simulations 
#                        are taking place over multiple nodes so that file names 
#                        should include the name of the current node
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

