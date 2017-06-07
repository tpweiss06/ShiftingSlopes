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
     
     # Define a matrix with the x and y movements for the eigth possible nearest
     #    neighbor movements according to this layout:
     #    1    2    3
     #    8    *    4
     #    7    6    5
     # The rows in the following matrix correspond to: [1,] numbered neighboring
     #    patch, [2,] required shift in x coordinates, [3,] required shift in y
     #    coordinates
     NearNeighbors <- matrix(NA, nrow = 3, ncol = 8)
     NearNeighbors[1,] <- 1:8
     NearNeighbors[2,] <- c(-1, 0, 1, 1, 1, 0, -1, -1)
     NearNeighbors[3,] <- c(1, 1, 1, 0, -1, -1, -1, 0)
     while(CurTime <= DispTime){
          # Update the current x and y coordinates for an individual
          patch <- sample(NearNeighbors[1,], size = 1)
          CurX[CurIndividual] <- CurX[CurIndividual] + NearNeighbors[2,patch]
          CurY[CurIndividual] <- CurY[CurIndividual] + NearNeighbors[3,patch]
          
          # Add in the next event time for the current individual and reset
          #    the current time and individual
          EventTimes[CurIndividual] <- rexp(n = 1, rate = 1) / DispTrait[CurIndividual,2]
          CurTime <- min(EventTimes)
          CurIndividual <- which(EventTimes == CurTime)
     }
     # Use the width of the landscape to correct any x coordinates outside
     #    of the allowed boundaries
     CurX <- ifelse( (CurX > width) | (CurX < 0), CurX %% width, CurX)
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

###### CalcEnvMean
# This function will use the relevant environmental parameters to calculate the
#    mean of the range capacity function over a given interval using the Mean
#    Value Theorem.
### INPUTS
# alpha, beta, gamma, and tau are all parameters in the range capacity function.
#    See the function write up for details on these parameters.
# a:      Lower interval bound
# b:      Upper interval bound
### OUTPUS
# This function will output a single value for the mean value of the function
#    over the given interval.
CalcEnvMean <- function(alpha, beta, gamma, tau, a, b){
     # First determine where the interval falls in the range, then calculate
     #    the appropriate mean (See write up)
     if(b <= beta){
          numerator <- exp(gamma * (b - beta + tau)) + 1
          denominator <- exp(gamma * (a - beta + tau)) + 1
          FullIntegral <- (alpha / gamma) * log(numerator / denominator)
          EnvMean <- (1 / (b - a)) * FullIntegral
     } else if(a >= beta){
          numerator <- exp(-1 * gamma * (b - beta - tau)) + 1
          denominator <- exp(-1 * gamma * (a - beta - tau)) + 1
          FullIntegral <- -1 * (alpha / gamma) * log(numerator / denominator)
          EnvMean <- (1 / (b - a)) * FullIntegral
     } else{
          PreNum <- exp(gamma * tau) + 1
          PreDen <- exp(gamma * (a - beta + tau)) + 1
          PostNum <- exp(-1 * gamma * (b - beta - tau)) + 1
          PostDen <- exp(gamma * tau) + 1
          FullIntegral <- (alpha / gamma) * log(PreNum / PreDen) + 
                         -1 * (alpha/gamma) * log(PostNum / PostDen)
          EnvMean <- (1 / (b - a)) * FullIntegral
     }
     return(EnvMean)
}

###### GetEnvQual
# This function will uses the CalcEnvMean function to extract a vector of 
#    environmental means for a given set of patches defined by their center
#    points and widths. 
### INPUTS
# alpha, beta, gamma, and tau are all parameters in the range capacity function.
#    See the function write up for details on these parameters.
# PatchCenters:     A vector of the center coordinates for discrete habitat
#                        patches
# PatchScale:       A single constant value used to modify the size of patches.
#                        By defaul this will be set to 1, but it could in 
#                        principle be used to explore the consequences of 
#                        discretizing space by letting it become arbitrarily
#                        small or alternatively it could allow a given range
#                        to contain less patches and thus save computational
#                        power or explore the effects of essentially lowering
#                        the population size.
### OUTPUTS
# This function will generate a vector of environmental quality values for each
#    patch as determined by the GetEnvMean function.
GetEnvQual <- function(alpha, beta, gamma, tau, PatchCenters, PatchScale = 1){
     
}


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

