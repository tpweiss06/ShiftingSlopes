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

###### rStudKern
# This function will perform a specified number of random draws from the one 
#    parameter student's t dispersal kernel used in the Shaw 2014 ppaer on the
#    consequences of risky dispersal.
### INPUTS
# n:      The number of random numbers to generate
# d:      The dispersal trait(s) for the individual(s) involved. This will be 
#              converted to another parameter to use in the probability 
#              distribution
### OUTPUTS
# n random numbers generated from the student's t dispersal kernel
rStudKern <- function(n, d){
     # First translate d to the u parameter used in the distribution
     u <- sqrt(d / 1.571)
     # Then generate n uniform random numbers between 0 and 1
     n_unif <- runif(n = n, min = 0, max = 1)
     # Now map them to the Student's t dispersal kernel using the CDF
     StudKernVals <- sqrt((u * n_unif) / (1 - n_unif))
     return(StudKernVals)
}

###### Disperse
# This function will use the previously calculated dispersal traits to calculate
#    new, post dispersal locations for individuals according to any of three
#    potential dispersal kernels: normal, exponential, or student's t dispersal
#    kernel (Shaw 2014).
### INPUTS:
# PopMat:      The population matrix to use
# DispTrait:   The two column matrix created by the CalcDispTrait function
# width:       The number of discrete cells defining the width of the word
#                   being simulated
# kern:        The dispersal kernel to be used. Can be "norm", "exp", or 
#                   "stud_t"
# PatchScale:  A constant value determining the width of individual patches on
#                   the Cartesian coordinate system defining both
#                   environmental quality and dispersal kernels.
### OUTPUTS:
# The function will return an updated population matrix with the post dispersal
#    columns filled in.
Disperse <- function(PopMat, DispTrait, width, kern, PatchScale = 1){
     # First determine the number of dispersal events to simulate
     N <- nrow(DispTrait)
     
     # Then generate the dispersal distances according to the type of kernel
     if(kern == "norm"){
          sigma <- (DispTrait[,2]^2 * pi) / 2
          dists <- abs(rnorm(n = N, mean = 0, sd = sigma))
     } else if(kern == "exp"){
          dists <- rexp(n = N, rate = 1 / DispTrait[,2])
     } else if(kern == "stud_t"){
          dists <- rStudKern(n = N, d = DispTrait[,2])
     }
     
     # Now generate the random, unbiased direction used by each individual
     angles <- runif(n = N, min = 0, max = 2 * pi)
     
     # Now calculate the new x and y coordinates of each individual as if their
     #    current location is the origin
     NewX <- dists * cos(angles)
     NewY <- dists * sin(angles)
     
     # Now use the PatchScale parameter to determine the number of x and y
     #    patches changed
     DeltaX <- floor(NewX / PatchScale)
     DeltaY <- floor(NewY / PatchScale)
     
     # Update the values in the population matrix
     PopMat[DispTrait[,1], "x1"] <- PopMat[DispTrait[,1], "x0"] + DeltaX
     PopMat[DispTrait[,1], "y1"] <- PopMat[DispTrait[,1], "y0"] + DeltaY
     
     # Finally, check and fix any x values that fall outside of the allowed
     #    width of the landscape and return the updated population matrix
     CurX <- PopMat[DispTrait[,1], "x1"]
     CurX <- ifelse( (CurX > width) | (CurX < 0), CurX %% width, CurX)
     CurX <- ifelse(CurX == 0, width, CurX)
     
     PopMat[DispTrait[,1], "x1"] <- CurX
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
# a:      Lower interval bound(s)
# b:      Upper interval bound(s)
### OUTPUS
# This function will output a single value for the mean value of the function
#    over the given interval.
CalcEnvMean <- function(alpha, beta, gamma, tau, a, b){
     # Determine the number of patches to calculate the mean for, then loop
     #    over them
     NumPatches <- length(a)
     EnvMean <- rep(NA, NumPatches)
     for(i in 1:NumPatches){
          # First determine where the interval falls in the range, then calculate
          #    the appropriate mean (See write up)
          if(b[i] <= beta){
               numerator <- exp(gamma * (b[i] - beta + tau)) + 1
               denominator <- exp(gamma * (a[i] - beta + tau)) + 1
               FullIntegral <- (alpha / gamma) * log(numerator / denominator)
               EnvMean[i] <- (1 / (b[i] - a[i])) * FullIntegral
          } else if(a[i] >= beta){
               numerator <- exp(-1 * gamma * (b[i] - beta - tau)) + 1
               denominator <- exp(-1 * gamma * (a[i] - beta - tau)) + 1
               FullIntegral <- -1 * (alpha / gamma) * log(numerator / denominator)
               EnvMean[i] <- (1 / (b[i] - a[i])) * FullIntegral
          } else{
               PreNum <- exp(gamma * tau) + 1
               PreDen <- exp(gamma * (a[i] - beta + tau)) + 1
               PostNum <- exp(-1 * gamma * (b[i] - beta - tau)) + 1
               PostDen <- exp(gamma * tau) + 1
               FullIntegral <- (alpha / gamma) * log(PreNum / PreDen) + 
                              -1 * (alpha/gamma) * log(PostNum / PostDen)
               EnvMean[i] <- (1 / (b[i] - a[i])) * FullIntegral
          }
     }
     return(EnvMean)
}

###### GetEnvQual
# This function uses the CalcEnvMean function to extract a vector of
#    environmental means for a given set of patches. 
### INPUTS
# alpha, beta, gamma, and tau are all parameters in the range capacity function.
#    See the function write up for details on these parameters.
# patches:          A vector of the patch numbers to get quality information for
# PatchScale:       A single constant value used to modify the size of patches.
#                        By default this will be set to 1, but it could in 
#                        principle be used to explore the consequences of 
#                        discretizing space by letting it become arbitrarily
#                        small or alternatively it could allow a given range
#                        to contain less patches and thus save computational
#                        power or explore the effects of essentially lowering
#                        the population size.
### OUTPUTS
# This function will generate a vector of environmental quality values for each
#    patch as determined by the GetEnvMean function.
GetEnvQual <- function(alpha, beta, gamma, tau, patches, PatchScale = 1){
     # Calculate the lower and upper bounds of each patch on the continuous scale
     #    by first calculating their center points
     centers <- patches * PatchScale
     lowers <- centers - PatchScale * 0.5
     uppers <- centers + PatchScale * 0.5
     
     # Now get and return the environmental quality score for each patch
     EnvQuals[i] <- CalcEnvMean(alpha = alpha, beta = beta, gamma = gamma, tau = tau, 
                                  a = lowers, b = uppers)
     return(EnvQuals)
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

