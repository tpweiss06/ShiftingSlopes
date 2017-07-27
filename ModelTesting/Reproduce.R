# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# Check the negative R error message
# Check appropriate lack of growth for very large population sizes
# Use width of 1 to allow for more precise control of Nt
# Check that setting R to 0 results in equal population size (on average)
# Check that setting K to 0 results in no growth (might have to add in a check to the main function for this)
# Check that moving the spatial location to more and more maladapted localities results in reduced growth on average

AllCols <- PopMatColNames(nFit = 1, nDisp = 1, monoecious = TRUE)
FitCols <- grep("^fit", AllCols)
DispCols <- grep("^disp", AllCols)
PopMat <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                     PopSize = 50, BetaInit = 0, width = 1, FitInit = 0.25,
                     FitDiv = 0, DispInit = 1, DispDiv = 1)
TraitMat <- CalcTraits(population = 1:50, PopMat = PopMat, FitColumns = FitCols,
                       DispColumns = DispCols)

# Set necessary parameters
LocalSel <- 2
tau <- 100
beta <- 0
omega <- 0.25
gamma <- 1
R0 <- 2
K0 <- 100
U <- c(0.25, 0.25)
Vm <- c(0.1, 0.1)
names(U) <- c("fit", "disp")
names(Vm) <- c("fit", "disp")
PatchScale <- 1

# Test that a useful error message is provided
Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = R0, K0 = K0,
          U = U, Vm = Vm, LocalSel = LocalSel, traits = TraitMat, PopMat = PopMat,
          EnvGradType = "k", ColumnNames = AllCols, FitColumns = FitCols, 
          DispColumns = DispCols, PatchScale = PatchScale)

Check1 <- TRUE

# Check the negative R error message
Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = -2, K0 = K0,
          U = U, Vm = Vm, LocalSel = LocalSel, traits = TraitMat, PopMat = PopMat,
          EnvGradType = "K", ColumnNames = AllCols, FitColumns = FitCols, 
          DispColumns = DispCols, PatchScale = PatchScale)

Check2 <- TRUE

# Check appropriate lack of growth for very large population sizes
LargePopSize <- ceiling(K0 * (1 + 1 / R0))
LargePopMat <- Initialize(ColumnNames = AllCols, FitCols = FitCols, DispCols = DispCols,
                     PopSize = LargePopSize, BetaInit = 0, width = 1, FitInit = 1,
                     FitDiv = 1, DispInit = 1, DispDiv = 1)
LargeTraitMat <- CalcTraits(population = 1:LargePopSize, PopMat = LargePopMat, FitColumns = FitCols,
                       DispColumns = DispCols)
# The reproduce function should return an empty matrix with 0 rows, signifying
#    an extinct population
Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = R0, K0 = K0,
          U = U, Vm = Vm, LocalSel = LocalSel, traits = LargeTraitMat, PopMat = LargePopMat,
          EnvGradType = "K", ColumnNames = AllCols, FitColumns = FitCols, 
          DispColumns = DispCols, PatchScale = PatchScale)

Check3 <- TRUE

# Check that setting R to 0 results in equal population size (on average) when
#    all individuals
TestPopSize <- rep(NA, 100)
for(i in 1:length(TestPopSize)){
     TestPop <- Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = 0, K0 = K0,
                          U = U, Vm = Vm, LocalSel = LocalSel, traits = TraitMat, PopMat = PopMat,
                          EnvGradType = "K", ColumnNames = AllCols, FitColumns = FitCols, 
                          DispColumns = DispCols, PatchScale = PatchScale)
     TestPopSize[i] <- nrow(TestPop)
}
# The following two quantities should be approximately equal
mean(TestPopSize) 
nrow(PopMat)

Check4 <- TRUE

# Check that setting K to 0 results in no growth
Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = R0, K0 = 0,
          U = U, Vm = Vm, LocalSel = LocalSel, traits = TraitMat, PopMat = PopMat,
          EnvGradType = "K", ColumnNames = AllCols, FitColumns = FitCols, 
          DispColumns = DispCols, PatchScale = PatchScale)

Check5 <- TRUE

# Check that moving the spatial location to more and more maladapted localities 
#    results in reduced growth on average
TestPopSize1 <- rep(NA, 100)
for(i in 1:length(TestPopSize1)){
     TestPop <- Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = R0, K0 = K0,
                          U = U, Vm = Vm, LocalSel = LocalSel, traits = TraitMat, PopMat = PopMat,
                          EnvGradType = "K", ColumnNames = AllCols, FitColumns = FitCols, 
                          DispColumns = DispCols, PatchScale = PatchScale)
     TestPopSize1[i] <- nrow(TestPop)
}


PopMat[,"x1"] <- rep(50, nrow(PopMat))
TestPopSize2 <- rep(NA, 100)
for(i in 1:length(TestPopSize2)){
     TestPop <- Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = R0, K0 = K0,
                          U = U, Vm = Vm, LocalSel = LocalSel, traits = TraitMat, PopMat = PopMat,
                          EnvGradType = "K", ColumnNames = AllCols, FitColumns = FitCols, 
                          DispColumns = DispCols, PatchScale = PatchScale)
     TestPopSize2[i] <- nrow(TestPop)
}

PopMat[,"x1"] <- rep(100, nrow(PopMat))
TestPopSize3 <- rep(NA, 100)
for(i in 1:length(TestPopSize3)){
     TestPop <- Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = R0, K0 = K0,
                          U = U, Vm = Vm, LocalSel = LocalSel, traits = TraitMat, PopMat = PopMat,
                          EnvGradType = "K", ColumnNames = AllCols, FitColumns = FitCols, 
                          DispColumns = DispCols, PatchScale = PatchScale)
     TestPopSize3[i] <- nrow(TestPop)
}

PopMat[,"x1"] <- rep(150, nrow(PopMat))
TestPopSize4 <- rep(NA, 100)
for(i in 1:length(TestPopSize4)){
     TestPop <- Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = R0, K0 = K0,
                          U = U, Vm = Vm, LocalSel = LocalSel, traits = TraitMat, PopMat = PopMat,
                          EnvGradType = "K", ColumnNames = AllCols, FitColumns = FitCols, 
                          DispColumns = DispCols, PatchScale = PatchScale)
     TestPopSize4[i] <- nrow(TestPop)
}

PopMat[,"x1"] <- rep(200, nrow(PopMat))
TestPopSize5 <- rep(NA, 100)
for(i in 1:length(TestPopSize5)){
     TestPop <- Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = R0, K0 = K0,
                          U = U, Vm = Vm, LocalSel = LocalSel, traits = TraitMat, PopMat = PopMat,
                          EnvGradType = "K", ColumnNames = AllCols, FitColumns = FitCols, 
                          DispColumns = DispCols, PatchScale = PatchScale)
     TestPopSize5[i] <- nrow(TestPop)
}

PopMat[,"x1"] <- rep(1000, nrow(PopMat))
TestPopSize6 <- rep(NA, 100)
for(i in 1:length(TestPopSize6)){
     TestPop <- Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega, R0 = R0, K0 = K0,
                          U = U, Vm = Vm, LocalSel = LocalSel, traits = TraitMat, PopMat = PopMat,
                          EnvGradType = "K", ColumnNames = AllCols, FitColumns = FitCols, 
                          DispColumns = DispCols, PatchScale = PatchScale)
     TestPopSize6[i] <- nrow(TestPop)
}

Check6 <- TRUE

if(Check1 & Check2 & Check3 & Check4 & Check5 & Check6){
     writeLines("Pass", con = "ModelTesting/Reproduce.txt")
}




