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

###### CalcTraits
# This function will calculate the fitness and dispersal trait values for a set 
#    of individuals based on their alleles
### INPUTS
# population:  A vector of the population (i.e. row numbers) for which to 
#                   calculate a trait values.
# PopMat:      The population matrix to use
# nFit:        The number of fitness loci in the current simulation
# nDisp:       The number of dispersal loci in the current simulation
### OUTPUTS
# The function will create a matrix with three columns. The first column will 
#    correspond to the ID (row number) of the individual, the second column
#    will be the fitness trait value, and the third column will be the dispersal
#    trait. If the traits are calculated for the entire population, as will 
#    normally be the case, then the first column will simply correspond to the 
#    row number, but including this column allows for the calculation of 
#    traits for a subset of the current population.
CalcTraits <- function(population, PopMat, FitColumns, DispColumns){
     # First create an empty matrix for the trait values and put
     #    the population IDs in the first column
     traits <- matrix(NA, nrow = length(population), ncol = 3)
     traits[,1] <- population
     
     # Now calculate the sum of each individual's quantitative loci for each 
     #    trait
     FitLociSum <- rowSums(PopMat[population, FitColumns])
     DispLociSum <- rowSums(PopMat[population, DispColumns])
     # Finally, store the trait values in the matrix and return it. 
     #    Note the exponential transformation of the dispersal trait to ensure
     #    that the trait value (mean distance dispersed) is bounded below by 0.
     traits[,2] <- FitLociSum
     traits[,3] <- exp(DispLociSum)
     colnames(traits) <- c("ID", "fit", "disp")
     return(traits)
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
# traits:      The three column matrix created by the CalcTrait function
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
Disperse <- function(PopMat, traits, width, kern, PatchScale = 1){
     # First determine the number of dispersal events to simulate
     N <- nrow(traits)
     
     # Then generate the dispersal distances according to the type of kernel
     if(kern == "norm"){
          sigma <- (traits[,"disp"]^2 * pi) / 2
          dists <- abs(rnorm(n = N, mean = 0, sd = sigma))
     } else if(kern == "exp"){
          dists <- rexp(n = N, rate = 1 / traits[,"disp"])
     } else if(kern == "stud_t"){
          dists <- rStudKern(n = N, d = traits[,"disp"])
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
     
     # Finally, check and fix any y values that fall outside of the allowed
     #    width of the landscape and return the updated population matrix
     CurY <- PopMat[traits[,"ID"], "y1"]
     CurY <- ifelse( (CurY > width) | (CurY < 0), CurY %% width, CurY)
     CurY <- ifelse(CurY == 0, width, CurY)
     
     PopMat[traits[,"ID"], "y1"] <- CurY
     return(PopMat)
}

###### FullSim
# Currently this is only mean to stand in as a template so that I can think
#    through the rest of the functions I will need.
### INPUTS
### OUTPUTS
FullSim <- function(){
     # First, create the column names for the population matrices used in this
     #    simulation and store the indices for the fitness and dispersal columns
     ColumnNames <- PopMatColNames(nFit = nFit, nDisp = nDisp, 
                                   monoecious = monoecious)
     FitCols <- grep("^fit", ColumnNames)
     DispCols <- grep("^disp", ColumnNames)
     
     # Next initialize the generation 0 founding population
     PopMat <- initialize(ColumnNames = ColumnNames, FitCols = FitCols, 
                          DispCols = DispCols)
     
     # Set up objects to hold whatever summary statistics we decide on here and
     #    populate them with the initial population values
     
     # Calculate the time points for the end of the shifting and the total time
     EndShift <- BurnIn + LengthShift
     TotalTime <- BurnIn + LengthShift + BurnOut
     
     # Calculate the beta values for during the period of climate change
     BetaShift <- ChangeClimate(BetaInit, LengthShift, ClimSpeed)
     # Now run through the actual simulation
     for(g in 1:TotalTime){
          if(g <= BurnIn){
               beta <- BetaInit
          } else if( (g > BurnIn) & (g <= EndShift) ){
               beta <- BetaShift[g - BurnIn]
          } else if( (g > EndShift) & (g <= TotalTime) ){
               beta <- BetaShift[LengthShift]
          }
          RepPopMat <- Reproduce(PopMat, beta = beta)
          DispPopMat <- Disperse(RepPopMat)
          PopMat <- DispPopMat
          
          # Keep track of all summary statistics here
     }
     
     # Finally, save the results here
     CurDirectory <- getwd()
     ResultsDir <- GetSafeID(ParentDirectory = CurDirectory, parallel = FALSE)
     dir.create(ResultsDir)
     write.csv(PopMat, file = paste(ResultsDir, "PopMat.csv", sep = "/"), 
               row.names = FALSE, quote = FALSE)
     SaveParams(FilePath = paste(ResultsDir, "Params.R", sep = "/"))
     # Either save or return my summary statistics here
     return(NULL)
     # I need to think carefully about a useful thing for the function to return
     #    (or if there is one).
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

###### initialize
# This function will initialize a population matrix of founders to start a 
#    simulation.
### INPUTS
# ColumnNames: The column names for the population matrix in the current 
#                   simulation
# FitCols:     The column indices for the fitness loci
# DispCols:    The column indices for the dispersal loci
# PopSize:     The number of founders to create in the intitial population
# BetaInit:    The starting location of the range center where individuals will
#                   initialize
# SexRatio:    The sex ratio of the founding population (if dioecious). Set to
#                   0.5 by default
# FitInit:     The mean initial value for all fitness loci
# FitDiv:      The standard deviation of the initial distribution of fitness 
#                   loci
# DispInit:    The mean initial value for all dispersal loci
# DispDiv:     The standard deviation of the initial distribution of dispersal 
#                   loci
### OUTPUS
# A filled in population matrix to start generation 0
initialize <- function(ColumnNames, FitCols, DispCols, PopSize, BetaInit, 
                       SexRatio = 0.5, FitInit, FitDiv, DispInit, DispDiv, ...){
     # First make an empty population matrix with the correct names
     PopMat <- matrix(NA, nrow = PopSize, ncol = length(ColumnNames))
     colnames(PopMat) <- ColumnNames
     
     # Next, fill in the x1 and y1 columns for the founders (the founders are
     #    considered post dispersal and will reproduce next)
     PopMat[,"x1"] <- BetaInit
     PopMat[,"y1"] <- sample(1:width, size = PopSize, replace = TRUE)
     
     # Fill in the sex column if it is present
     if("sex" %in% ColumnNames){
          PopMat[,"sex"] <- rbinom(n = PopSize, size = 1, prob = SexRatio)
     }
     
     # Now fill in the fitness and dispersal allele columns
     PopMat[,FitCols] <- rnorm(n = PopSize * length(FitCols), mean = FitInit,
                               sd = FitDiv)
     PopMat[,DispCols] <- rnorm(n = PopSize * length(DispCols), mean = DispInit,
                                sd = DispDiv)
     
     # Return the filled in initial population matrix
     return(PopMat)
}

###### PopMatColNames
# This function generates the column names to use for a given simulation
### INPUTS
# nFit:        the number of loci defining an individual's fitness
# nDisp:       the number of loci defining an individual's dispersal ability
# monoecious:  a boolean value indicating whether the individuals in the simulation
#              are monoecious (i.e. possessing both the female and male 
#              reproductive organs in the same individual) or not.
# example:     a boolean value indicating whether the function should simply 
#                   return an example character vector containing the column 
#                   names and order used
### OUTPUTS
# This function creates a single, empty matrix wth PopSize rows and nFit + nDisp
#    + 5 columns with appropriate names (see function definition for details). 
#    If example = TRUE, however, the function simply returns a simple character
#    vector with the column names and layout.
PopMatColNames <- function(nFit, nDisp, monoecious, example = FALSE){
     if(example){
          MatNames <- c("x0", "y0", "x1", "y1", "sex", "fit1", "...", "fitN",
                        "disp1", "...", "dispN")
          return(MatNames)
     }
     # First determine the number of columns and make the population matrix
     MatCols <- ifelse(monoecious, (nFit + nDisp + 4), (nFit + nDisp + 5) )
     
     # Next create an empty vector of column names
     MatNames <- rep(NA, MatCols)
     
     # Create the names of the first 4 or 5 columns depending on whether
     #    individuals are monoecious and therefore whether a sex column is
     #    needed.
     if(monoecious){
          TraitStart <- 4
          MatNames[1:TraitStart] <- c("x0", "y0", "x1", "y1")
     } else{
          TraitStart <- 5
          MatNames[1:TraitStart] <- c("x0", "y0", "x1", "y1", "sex")
     }
     
     # Next create the column names for the fitness and dispersal traits
     for(n in 1:nFit){
          ColName <- paste("fit", n, sep = "")
          MatNames[TraitStart + n] <- ColName
     }
     for(n in 1:nDisp){
          ColName <- paste("disp", n, sep = "")
          MatNames[TraitStart + nFit + n] <- ColName
     }
     
     # Return the names for the population matrix
     return(MatNames)
}

###### GetSafeID
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
GetSafeID <- function(ParentDirectory, parallel = FALSE){
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

