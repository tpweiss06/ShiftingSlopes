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
     PopMat[traits[,"ID"], "x1"] <- PopMat[traits[,"ID"], "x0"] + DeltaX
     PopMat[traits[,"ID"], "y1"] <- PopMat[traits[,"ID"], "y0"] + DeltaY
     
     # Finally, check and fix any y values that fall outside of the allowed
     #    width of the landscape and return the updated population matrix
     CurY <- PopMat[traits[,"ID"], "y1"]
     CurY <- ifelse( (CurY > width) | (CurY < 0), CurY %% width, CurY)
     CurY <- ifelse(CurY == 0, width, CurY)
     
     PopMat[traits[,"ID"], "y1"] <- CurY
     return(PopMat)
}

###### RelFit
# This function will calculate the relative fitness values for individuals 
#    depending on their position within the range and their respective phenotype
#    values.
### INPUTS
# LocalSel:    A value determining the strength of local selection across the
#                   range. Higher values will correspond to a steeper gradient
#                   in optimum phenotype values from one end of the range to the
#                   other.
# tau:         The tau value used to define the width of the range
# beta:        The beta value used to define the center of the range
# traits:      The trait matrix calculated by the CalcTraits function
# PopMat:      The population matrix containing the focal individuals
# individuals: A vector of the individual IDs (row numbers) of the focal 
#                   individuals.
# omega:       This parameter is the inverse of the strength of stabilizing 
#                   selection.
### OUTPUTS
# This function will generate a vector of relative fitness values for each
#    individual in the same order they are presented in the phenotype vector.
RelFit <- function(LocalSel, tau, beta, traits, PopMat, individuals, omega){
     # Calculate the "slope" of the phenotypic optimum function using tau
     #    and the LocalSel argument
     lambda <- LocalSel / tau
     
     # Calculate the Zopt values for each individual depending on where they
     #    are in space.
     numerator <- exp(lambda * (PopMat[individuals, "x1"] - beta))
     denominator <- 1 + exp(lambda * (PopMat[individuals, "x1"] - beta))
     Zopt <- numerator / denominator
     
     # Finally, use the equation for stabilizing selection to calculate each
     #    individual's relative fitness and return that as a vector
     RelFits <- exp(-1 * (traits[individuals, "fit"] - Zopt)^2 / (2*omega^2))
     return(RelFits)
}

###### MatFill
# This function will fill in the appropriate values for a new population matrix
#    in the Reproduce function.
### INPUTS
# RealizedNtp1:      A vector of the stochastically determined population size for
#                        each occupied patch in the next generation
# ColumnNames:       A vector of the column names being used in the population
#                        matrix for these simulations
# OccPatches:        A two column matrix with the x and y coordinates of all
#                        occupied patches in the current generation.
# SexRatio:          The sex ratio of the population (set to 0.5 by default)
# RelFits:           A vector of the relative fitness values of individuals in the
#                        population.
# U:                 The diploid mutation rate for quantitative trait z
# Vm:                The mutational variance component of the total quantitative
#                        genetic variance
# FitColumns:        A vector of the columns corresponding to fitness traits
# DispColumns:       A vector of the columns corresponding to dispersal
# PopMat:            The population matrix for the current generation
### OUTPUTS
# This function will return a new matrix composed of all the offspring and their
#    location, loci, etc. for the next generation.
MatFill <- function(RealizedNtp1, ColumnNames, OccPatches, SexRatio, RelFits,
                    U, Vm, FitColumns, DispColumns, PopMat){
     # Create the new matrix for the next population
     NewPopMat <- matrix(NA, nrow = sum(RealizedNtp1), ncol = length(ColumnNames))
     colnames(NewPopMat) <- ColumnNames
     
     # Filter for only the patches that produced offspring
     NewOccPatches <- which(RealizedNtp1 != 0)
     Ntp1 <- RealizedNtp1[NewOccPatches]
     
     # Now step through and fill in the matrix
     for(i in 1:length(NewOccPatches)){
          # Calculate the starting and ending row for offspring in the current
          #    patch
          StartRow <- ifelse(i == 1, 1, 1 + sum(Ntp1[1:(i-1)]))
          EndRow <- sum(Ntp1[1:i])
          
          # Fill in the location details
          NewPopMat[StartRow:EndRow, "x0"] <- OccPatches[NewOccPatches[i], "x1"]
          NewPopMat[StartRow:EndRow, "y0"] <- OccPatches[NewOccPatches[i], "y1"]
          
          # Depending on life history being modeled (monoecious vs. dioecious),
          #    add in sex information if relevant and determine parents for each
          #    offspring
          if("sex" %in% ColumnNames){
               NewPopMat[StartRow:EndRow, "sex"] <- rbinom(n = Ntp1[i], size = 1,
                                                       prob = SexRatio)
               
               # Identify the females and males present in the current patch
               females <- which( (PopMat[,"sex"] == 1) & 
                              (PopMat[,"x1"] == OccPatches[NewOccPatches[i], "x1"]) &
                              (PopMat[,"y1"] == OccPatches[NewOccPatches[i], "y1"]))
               males <- which( (PopMat[,"sex"] == 0) & 
                              (PopMat[,"x1"] == OccPatches[NewOccPatches[i], "x1"]) &
                              (PopMat[,"y1"] == OccPatches[NewOccPatches[i], "y1"]))
               
               # Now extract the relative fitness values for the males and
               #    females
               FemaleFits <- RelFits[females]
               MaleFits <- RelFits[males]
                    
               # Sample from the pool of females and males according to their
               #    relative fitnesses to select parents for each offspring, but
               #    only use the sample command if there are more than one male
               #    or female available to avoid undesired behavior from the
               #    sample function
               if(length(females) == 1){
                    parent1 <- rep(females, Ntp1[i])
               } else{
                    parent1 <- sample(females, size = Ntp1[i], replace = TRUE,
                                      prob = FemaleFits)
               }
               if(length(males) == 1){
                    parent2 <- rep(males, Ntp1[i])
               } else{
                    parent2 <- sample(males, size = Ntp1[i], replace = TRUE,
                                      prob = MaleFits)
               }
               parents <- cbind(parent1, parent2)
          } else{
               # Define the pool of potential parents for the offspring and 
               #    their associated relative fitnesses
               ParentPool <- which( (PopMat[,"x1"] == OccPatches[NewOccPatches[i], "x1"]) &
                              (PopMat[,"y1"] == OccPatches[NewOccPatches[i], "y1"]))
               ParentFits <- RelFits[ParentPool]
               
               # Next select both parents, again avoiding unwanted behavior
               #    from the sample function
               if(length(ParentPool) == 1){
                    parent1 <- rep(ParentPool, Ntp1[i])
                    parent2 <- rep(ParentPool, Ntp1[i])
               } else if(length(ParentPool) == 2){
                    parent1 <- rep(ParentPool[1], Ntp1[i])
                    parent2 <- rep(ParentPool[2], Ntp1[i])
               } else{
                    parent1 <- sample(ParentPool, size = Ntp1[i], replace = TRUE, 
                                      prob = ParentFits)
                    parent2 <- rep(NA, Ntp1[i])
                    for(j in 1:Ntp1[i]){
                         NewParentPool <- setdiff(ParentPool, parent1[j])
                         parent2[j] <- sample(NewParentPool, size = 1,
                                              prob = RelFits[NewParentPool])
                    }
               }
               parents <- cbind(parent1, parent2)
          }
          # Fill in the loci according to the inheritence function
          NewPopMat[StartRow:EndRow, FitColumns] <- Inheritence(Cols = FitColumns, 
                                             parents = parents, PopMat = PopMat,
                                             U = U["fit"], Vm = Vm["fit"])
          NewPopMat[StartRow:EndRow, DispColumns] <- Inheritence(Cols = DispColumns, 
                                             parents = parents, PopMat = PopMat,
                                             U = U["disp"], Vm = Vm["disp"])
     }
     return(NewPopMat)
}

###### Inheritence
# This function performs inheritence of alleles from both parents, assuming
#    independent segregation and no crossover events. It incorporates mutations
#    in alleles via the U and Vm parameters (see below).
# INPUTS
# parents:     A matrix with two columns and a number of rows equal to the 
#                   number of offpsring being created. Each row contains the
#                   IDs (row indices in PopMat) of the parent(s) of the 
#                   offspring corresponding to that row.
# Cols:        A vector containing the focal columns containing the loci for the
#                   trait being inherited.
# PopMat:      The matrix containing all population information for the parental
#                   generation.
# U:           The diploid mutation rate for the focal quantitative trait
# Vm:          The mutational variance component of the total quantitative
#                   genetic variance
# OUTPUTS
# This function will return a matrix consisting of a number of rows equal to
#    the number of offspring produced and a number of columns equal to the 
#    number of loci defining the trait under consideration.
Inheritence <- function(Cols, parents, PopMat, U, Vm){
     NumOffspring <- nrow(parents)
     NumLoci <- length(Cols) / 2
     SegregatedLoci <- matrix(NA, nrow = NumOffspring, ncol = 2*NumLoci)
     for(i in 1:nrow(parents)){
          ParentLoci <- PopMat[parents[i,], Cols]
          Parent1Alleles <- 1:NumLoci + NumLoci * sample(c(0,1), replace = TRUE, size = NumLoci)
          Parent2Alleles <- 1:NumLoci + NumLoci * sample(c(0,1), replace = TRUE, size = NumLoci)
          SegregatedLoci[i,] <- c(ParentLoci[1,Parent1Alleles], ParentLoci[2,Parent2Alleles])
     }
     
     # Calculate the per locus mutation probability and the standard deviation
     #    for the mutational effects. The following calculations are taken from
     #    Gilbert et al. 2017 AmNat
     PerLocusProb <- U / (2 * NumLoci)
     MutStd <- sqrt(Vm / U)
     
     # Use the per locus mutation probability to identify a number of mutations
     #    for each new offspring
     NumMut <- rbinom(n = NumOffspring, size = NumLoci, prob = PerLocusProb)
     
     # Now step through each offspring in which a mutation takes place and alter
     #    allele values appropriately
     MutOffspring <- which(NumMut != 0)
     for(i in MutOffspring){
          MutLocus <- sample(1:NumLoci, size = NumMut[i])
          SegregatedLoci[i,MutLocus] <- rnorm(mean = SegregatedLoci[i,MutLocus], 
                                              sd = MutStd, n = NumMut[i])
     }
     
     # return the segregated and mutated loci
     return(SegregatedLoci)
}

###### Reproduce
# This function will use the previously calculated fitness traits and perform
#    logisitic population growth in all occupied patches with offspring assigned
#    to parental pairs according to individual fitness values.
### INPUTS:
# alpha, beta, gamma, tau:    Parameters corresponding to the range capacity fxn
# omega:            The inverse of the strength of stabilizing selection.             
# R0, K0:           The growth and carrying capacity parameters for logistic
#                        growth. These values are the maximum attainable values
#                        which are then used to compute realized values
#                        throughout the range.
# U:                A vector of 2 entries with the diploid mutation rate for
#                        each trait.
# Vm:               A vector of 2 entries with the mutational variance for each
#                        trait.
# LocalSel:         A value determining the gradient in local phenotypic optima
#                        across the range.
# traits:           A two column matrix of trait values created by CalcTraits
# PopMat:           The matrix with all information on the current generation
# EnvGradType:      The type of environmental gradient being measured (either "K"
#                        or "R")
# ColumnNames:      A vector of the column names relevant for the simulation.
# SexRatio:         The sex ratio of the population if individuals are dioecious.
#                        By default, this is set to 0.5
# FitColumns:       A vector of the column indices corresponding to fitness loci
# DispColumns:      A vector of the column indices corresponding to dispersal loci
### OUTPUTS:
# The function will return new population matrix with the offspring resulting
#    from reproduction with their x0 and y0 columns corresponding to their
#    natal patches.
Reproduce <- function(alpha = 1, beta, gamma, tau, omega, R0, K0, U, Vm, LocalSel, 
                      traits, PopMat, EnvGradType, ColumnNames, SexRatio = 0.5,
                      FitColumns, DispColumns, PatchScale){
     # Check that R0 is non-negative before anything else
     if(R0 < 0){
          write("Negative values are not recomended for R0 as they can result in undesired behavior in the logistic growth equation",
                stderr())
          return(NULL)
     }
     
     # First calculate the relative fitness of each individual in the population
     RelFits <- RelFit(LocalSel = LocalSel, tau = tau, beta = beta, traits = traits, 
                       PopMat = PopMat, individuals = 1:nrow(PopMat), omega = omega)
     
     # Next generate a list of all the occupied patches and calculate the
     #    environmental quality score for each occupied patch
     OccPatches <- unique(PopMat[,c("x1","y1")])
     PatchEnvQual <- GetEnvQual(alpha = alpha, beta = beta, gamma = gamma, 
                              tau = tau, patches = OccPatches[,"x1"], 
                              PatchScale = PatchScale)
     
     # Create objects to hold the current population and the expected size for
     # the next generation's population size
     Ntp1 <- rep(NA, nrow(OccPatches))
     Nt <- vector(mode = "list", length = nrow(OccPatches))
     
     # Loop through each occupied patch and calculate the expected population
     #    size in the next generation
     for(i in 1:nrow(OccPatches)){
          # Extract the quality of the current patch and the population size
          PatchQual <- PatchEnvQual[i]
          Nt[[i]] <- which( (PopMat[,"x1"] == OccPatches[i,"x1"]) & 
                          (PopMat[,"y1"] == OccPatches[i,"y1"]) )
          PopSize <- length(Nt[[i]])
          
          # Calculate the relevant demographic parameter and return a helpful
          #    error message if needed
          if(EnvGradType == "K"){
               PatchK <- PatchQual * K0 * PatchScale
               PatchR <- R0
          } else if(EnvGradType == "R"){
               PatchR <- PatchQual * R0
               PatchK <- K0 * PatchScale
          } else{
               write("Invalid type of environmental gradient", stderr())
               return(NULL)
          }
          
          # Calculate the expected population size based on the logistic equation,
          #    correcting for cases where a large ratio of Nt to K will mathematically
          #    lead to unrealistically large increases in population size. Also
          #    correcting for the times when both R and K are 0, leading to a 
          #    NaN result.
          if( (PatchR == 0) & (PatchK == 0)){
               ExpPopGrowth <- 0
          } else if(PopSize < (PatchK * (1 + 1/PatchR))){
               ExpPopGrowth <- 1 + PatchR * (1 - PopSize / PatchK)
          } else{
               ExpPopGrowth <- 0
          }
          
          # Now calculate the mean fitness of the current patch and use it to
          #    calculate the expected population growth after accounting for 
          #    hard selection on local fitness
          CurPop <- which((PopMat[,"x1"] == OccPatches[i, "x1"]) &
                               (PopMat[,"y1"] == OccPatches[i,"y1"]))
          MeanFit <- mean(RelFits[CurPop])
          AdjPopGrowth <- MeanFit * ExpPopGrowth * PopSize
          
          # If using dioecious individuals, ensure there are at least 1 male and
          #    female to allow for reproduction (mate finding Allee effect).
          if("sex" %in% ColumnNames){
               NumFemales <- length(which( (PopMat[,"sex"] == 1) & 
                                   (PopMat[,"x1"] == OccPatches[i, "x1"]) &
                                   (PopMat[,"y1"] == OccPatches[i, "y1"])))
               NumMales <- length(which( (PopMat[,"sex"] == 0) & 
                                   (PopMat[,"x1"] == OccPatches[i, "x1"]) &
                                   (PopMat[,"y1"] == OccPatches[i, "y1"])))
               if( (NumFemales >= 1) & (NumMales >= 1) ){
                    Ntp1[i] <- AdjPopGrowth
               } else{
                    Ntp1[i] <- 0
               }
          } else{
               Ntp1[i] <- AdjPopGrowth
          }
     }
     
     # Check for negative values of Ntp1 and reset those to 0
     LambdaVals <- ifelse(Ntp1 < 0, 0, Ntp1)
     
     # Use the expected population sizes to generate the realized population 
     #    sizes for each patch
     RealizedNtp1 <- rpois(n = length(Ntp1), lambda = LambdaVals)
     
     # Now make and return a new population matrix for the offspring, if there
     #    are any
     if(sum(RealizedNtp1) > 0){
          Ntp1Mat <- MatFill(RealizedNtp1 = RealizedNtp1, ColumnNames = ColumnNames,
                             OccPatches = OccPatches, SexRatio = SexRatio, PopMat = PopMat,
                             RelFits = RelFits, FitColumns = FitColumns, 
                             DispColumns = DispColumns, U = U, Vm = Vm)
     } else{
          Ntp1Mat <- matrix(NA, nrow = 0, ncol = length(ColumnNames))
          colnames(Ntp1Mat) <- ColumnNames
     }
     return(Ntp1Mat)
}

# ------------------------------------------------------------------------------
# -------------------------- Environmental Functions ---------------------------
# ------------------------------------------------------------------------------

###### ChangeClimate
# This function will create a vector of beta values depending on how fast the 
#    climate is set to change, assuming a constant rate of change
### INPUTS
# BetaInit:         The starting value of beta
# LengthShift:      The length of the climate shift
# ClimSpeed:        The speed with which the climate shifts
### OUTPUT
# A one dimensional vector containing the shifted beta values at each time point
#    for the duration of the climate change period
ChangeClimate <- function(BetaInit, LengthShift, ClimSpeed){
     # First check for a negative speed
     if(ClimSpeed < 0){
          write("Negative speed is not supported for climate change", stderr())
          return(NULL)
     }
     # Next check for a 0 value for the length of climate change
     if(LengthShift == 0){
          return(NULL)
     }
     # Perform the climate shift for the range center
     TimeVec <- 1:LengthShift
     BetaVec <- BetaInit + ClimSpeed*TimeVec
     return(BetaVec)
}

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
     EnvQuals <- CalcEnvMean(alpha = alpha, beta = beta, gamma = gamma, tau = tau, 
                                  a = lowers, b = uppers)
     return(EnvQuals)
}


# ------------------------------------------------------------------------------
# --------------------------- Bookkeeping Functions ----------------------------
# ------------------------------------------------------------------------------

###### Initialize
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
Initialize <- function(ColumnNames, FitCols, DispCols, PopSize, BetaInit, width, 
                       SexRatio = 0.5, FitInit, FitDiv, DispInit, DispDiv){
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
          MatNames <- c("x0", "y0", "x1", "y1", "sex", "fit1_1", "fit1_2", "...", "fit1_N",
                        "fit2_1", "...", "fit2_N", "disp1_1", "...", "disp2_N")
          return(MatNames)
     }
     # First determine the number of columns and make the population matrix
     MatCols <- ifelse(monoecious, (2*nFit + 2*nDisp + 4), (2*nFit + 2*nDisp + 5) )
     
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
          ColName <- paste("fit1", n, sep = "_")
          MatNames[TraitStart + n] <- ColName
     }
     for(n in 1:nFit){
          ColName <- paste("fit2", n, sep = "_")
          MatNames[TraitStart + nFit + n] <- ColName
     }
     for(n in 1:nDisp){
          ColName <- paste("disp1", n, sep = "_")
          MatNames[TraitStart + 2*nFit + n] <- ColName
     }
     for(n in 1:nDisp){
          ColName <- paste("disp2", n, sep = "_")
          MatNames[TraitStart + 2*nFit + nDisp + n] <- ColName
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

###### SaveParams
# This function will unpack a list of parameters and save them in an output file
#    that can be sourced to load all parameters into the environment. In this 
#    way, the parameters for each simulation are saved and easily accessible and
#    it allows the parameters for the simulations to be passed in to the main 
#    function in an easier format (i.e. a single list).
### INPUTS
# parameters:  a list containing all the parameter values necessary for a 
#                   simulation.
# FilePath:    a character vector for the results directory for the current
#                   simulation.
### OUTPUTS
# This function will unpack the parameter vector and create a .R file which
#    creates objects with the parameter names and values and named for the
#    specific simulation.
SaveParams <- function(parameters, FilePath){
     # First cheack that all necessary parameters are included
     ParamCheck <- names(parameters) == c("BetaInit", "gamma", "tau", "LocalSel", "omega",
                              "U", "Vm", "nFit", "nDisp", "R0", "K0", "width",
                              "kern", "EnvGradType", "monoecious", "BurnIn",
                              "BurnOut", "LengthShift", "ClimSpeed", "InitPopSize",
                              "FitInit", "FitDiv", "DispInit", "DispDiv", "PatchScale")
     if(sum(ParamCheck) != length(ParamCheck)){
          write("Incorrect names or number of input parameters", stderr())
          return(NULL)
     }
     
     OutFile <- paste(FilePath, "parameters.R", sep = "/")
     sink(OutFile)
     cat("BetaInit <- ", parameters$BetaInit, "\n", sep = "")
     cat("gamma <- ", parameters$gamma, "\n", sep = "")
     cat("tau <- ", parameters$tau, "\n", sep = "")
     cat("LocalSel <- ", parameters$LocalSel, "\n", sep = "")
     cat("omega <- ", parameters$omega, "\n", sep = "")
     cat("U <- c(", parameters$U[1], ",", parameters$U[2], ")\n", sep = "")
     cat("Vm <- c(", parameters$Vm[1], ",", parameters$Vm[2], ")\n", sep = "")
     cat("names(U) <- names(Vm) <- c(\"fit\", \"disp\")\n")
     cat("nFit <- ", parameters$nFit, "\n", sep = "")
     cat("nDisp <- ", parameters$nDisp, "\n", sep = "")
     cat("R0 <- ", parameters$R0, "\n", sep = "")
     cat("K0 <- ", parameters$K0, "\n", sep = "")
     cat("width <- ", parameters$width, "\n", sep = "")
     cat("kern <- \"", parameters$kern, "\"\n", sep = "")
     cat("EnvGradType <- \"", parameters$EnvGradType, "\"\n", sep = "")
     cat("monoecious <- ", parameters$monoecious, "\n", sep = "")
     cat("BurnIn <- ", parameters$BurnIn, "\n", sep = "")
     cat("BurnOut <- ", parameters$BurnOut, "\n", sep = "")
     cat("LengthShift <- ", parameters$LengthShift, "\n", sep = "")
     cat("ClimSpeed <- ", parameters$ClimSpeed, "\n", sep = "")
     cat("InitPopSize <- ", parameters$InitPopSize, "\n", sep = "")
     cat("FitInit <- ", parameters$FitInit, "\n", sep = "")
     cat("FitDiv <- ", parameters$FitDiv, "\n", sep = "")
     cat("DispInit <- ", parameters$DispInit, "\n", sep = "")
     cat("DispDiv <- ", parameters$DispDiv, "\n", sep = "")
     cat("PatchScale <- ", parameters$PatchScale,  "\n", sep = "")
     sink()
}


# ------------------------------------------------------------------------------
# ------------------------ Full Simulation Function ----------------------------
# ------------------------------------------------------------------------------

###### FullSim
# Currently this is only mean to stand in as a template so that I can think
#    through the rest of the functions I will need.
### INPUTS
# parameters:  A list with all necessary parameter values for a single model
#                   run. These will be unpacked and stored by the SaveParams
#                   function. It will include all the parameters necessary for
#                   the following functions as well as a few additional time
#                   keeping parameters: 
#                   Functions --
#                   PopMatColNames(), Initialize(), ChangeClimate(), Reproduce(), 
#                        Disperse(), and GetSafeID()
#                   Time keeping parameters --
#                   BurnIn:        Length of time pre shift
#                   LengthShift:   Duration of the shift
#                   BurnOut:       Length of time post shift
# parallel: A boolean variable indicating whether the simulations are being run
#              on a server or not (which affects how file paths are determined).
### OUTPUTS
FullSim <- function(parameters, parallel = FALSE){
     # First generate a save directory name and create it to save all output
     #    from the simulation
     CurDirectory <- getwd()
     ResultsDir <- GetSafeID(ParentDirectory = CurDirectory, parallel = parallel)
     dir.create(ResultsDir)
     
     # Next, save the parameters used for this simulation and source the file
     #    to have access to them within the function
     SaveParams(parameters = parameters, FilePath = ResultsDir)
     source(paste(ResultsDir, "parameters.R", sep = "/"))
     
     # Next, create the column names for the population matrices used in this
     #    simulation and store the indices for the fitness and dispersal columns
     ColumnNames <- PopMatColNames(nFit = nFit, nDisp = nDisp, 
                                   monoecious = monoecious)
     FitCols <- grep("^fit", ColumnNames)
     DispCols <- grep("^disp", ColumnNames)
     
     # Next initialize the generation 0 founding population and allow it to
     #    reproduce
     PopMat <- Initialize(ColumnNames = ColumnNames, FitCols = FitCols, 
                          DispCols = DispCols, PopSize = InitPopSize, 
                          BetaInit = BetaInit, FitInit = FitInit, FitDiv = FitDiv, 
                          DispInit = DispInit, DispDiv = DispDiv, width = width)
     traits <- CalcTraits(population = 1:InitPopSize, PopMat = PopMat, FitColumns = FitCols,
                          DispColumns = DispCols)
     PopMat <- Reproduce(beta = BetaInit, gamma = gamma, tau = tau, omega = omega,
                         R0 = R0, K0 = K0, U = U, Vm = Vm, LocalSel = LocalSel,
                         traits = traits, PopMat = PopMat, EnvGradType = EnvGradType,
                         ColumnNames = ColumnNames, FitColumns = FitCols, 
                         DispColumns = DispCols, PatchScale = PatchScale)
     
     # Calculate the time points for the end of the shifting and the total time
     EndShift <- BurnIn + LengthShift
     TotalTime <- BurnIn + LengthShift + BurnOut
     
     # Set up objects to hold summary information including patch abundance,
     #    patch genetic diversity in each trait, mean and standard deviations
     #    for patch trait values. As a first step, calculate number of patches
     #    to track for each of these objects which will be centered on the
     #    range center (beta)
     RangeK <- K0 * PatchScale * GetEnvQual(alpha = 1, beta = BetaInit, gamma = gamma, 
                         tau = tau, patches = -1000:1000, PatchScale = PatchScale)
     HalfRange <- sum(RangeK >= 1) %/% 2
     RangeLength <- 2 * (HalfRange + 10) + 1
     BetaIndex <- ceiling(RangeLength / 2)
     PatchAbund <- array(0, dim = c(TotalTime, width, RangeLength))
     PatchFitDiv <- array(NA, dim = c(TotalTime, width, RangeLength))
     PatchDispDiv <- array(NA, dim = c(TotalTime, width, RangeLength))
     PatchFitMean <- array(NA, dim = c(TotalTime, width, RangeLength))
     PatchDispMean <- array(NA, dim = c(TotalTime, width, RangeLength))
     
     # Calculate the beta values for during the period of climate change
     BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
                                ClimSpeed = ClimSpeed)
     # Now run through the actual simulation
     for(g in 1:TotalTime){
          if(g <= BurnIn){
               beta <- BetaInit
          } else if( (g > BurnIn) & (g <= EndShift) ){
               beta <- BetaShift[g - BurnIn]
          } else if( (g > EndShift) ){
               beta <- BetaShift[LengthShift]
          }
          
          # Check for extinction before dispersal and reproduction
          CurPopSize <- nrow(PopMat)
          if(CurPopSize != 0){
               traits <- CalcTraits(population = 1:CurPopSize, PopMat = PopMat, FitColumns = FitCols,
                                    DispColumns = DispCols)
               DispPopMat <- Disperse(PopMat = PopMat, traits = traits, width = width,
                                      kern = kern, PatchScale = PatchScale)
               RepPopMat <- Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega,
                                      R0 = R0, K0 = K0, U = U, Vm = Vm, LocalSel = LocalSel,
                                      traits = traits, PopMat = DispPopMat, EnvGradType = EnvGradType,
                                      ColumnNames = ColumnNames, FitColumns = FitCols, 
                                      DispColumns = DispCols, PatchScale = PatchScale)
               PopMat <- RepPopMat
               
               # Keep track of all summary statistics here
               OccPatches <- unique(PopMat[,c("x0","y0")])
               for(i in 1:nrow(OccPatches)){
                    CurIndices <- c(g, OccPatches[i,"y0"], (OccPatches[i,"x0"] - beta) +
                                         BetaIndex)
                    if( (CurIndices[3] < 1) | (CurIndices[3] > RangeLength) ){
                         write("Summary statistics not initialized with enough space", stderr())
                         return(NULL)
                    }
                    CurPop <- which((PopMat[,"x0"] == OccPatches[i, "x0"]) &
                                         (PopMat[,"y0"] == OccPatches[i,"y0"]))
                    PatchAbund[CurIndices] <- length(CurPop)
                    PatchFitDiv[CurIndices] <- sqrt(sd(PopMat[CurPop, FitCols]))
                    PatchFitMean[CurIndices] <- mean(PopMat[CurPop, FitCols])
                    PatchDispDiv[CurIndices] <- sqrt(sd(PopMat[CurPop, DispCols]))
                    PatchDispMean[CurIndices] <- mean(PopMat[CurPop, DispCols])
               }
          }
     }
     
     # Finally, save the results here
     write.csv(PopMat, file = paste(ResultsDir, "PopMat.csv", sep = "/"), 
               row.names = FALSE, quote = FALSE)
     save(PatchAbund, PatchFitDiv, PatchFitMean, PatchDispDiv, PatchDispMean,
          file = paste(ResultsDir, "SummaryStats.Rdata", sep = "/"))
     
     return(NULL)
}