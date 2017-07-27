# Set the working directory and load in the simulation functions
setwd("~/Desktop/RangeShifts/ShiftingSlopes/")
source("SimulationFunctions.R")

# To test the functions calculating the environmental conditions of the range, 
#    (1) first test the environmental means function by using a large sequence
#    of upper and lower bounds and examining the result. The function is symmetric
#    so the results should be as well and should generally resemble the shape of
#    the function. Then, (2) test the GenEnvQual function by giving it vectors
#    of unordered patch numbers and various scales for the patch width. Large
#    scales should be the average of smaller scales.

# Use the following parameters to check for symmetry
beta <- 0
tau <- 100
alpha <- 1
gamma <- 0.025
a <- seq(-150, 150, by = 10)
b <- seq(-140, 160, by = 10)

EnvMeans <- CalcEnvMean(alpha = alpha, beta = beta, gamma = gamma, tau = tau,
                             a = a, b = b)
plot(x = a, y = EnvMeans)

# Check that the patch scale argument works
SmallPatches <- GetEnvQual(alpha = alpha, beta = beta, gamma = gamma, tau = tau, 
                           patches = -150:150, PatchScale = 1)
MedPatches <- GetEnvQual(alpha = alpha, beta = beta, gamma = gamma, tau = tau, 
                         patches = -50:50, PatchScale = 3)
LargePatches <- GetEnvQual(alpha = alpha, beta = beta, gamma = gamma, tau = tau, 
                           patches = -15:15, PatchScale = 10)
plot(x = NA, y = NA, xlim = c(-150, 150), ylim = c(0, 1), xlab = "Spatial location",
     ylab = "Patch Quality", main = "")
points(x = seq(-150, 150, by = 10), y = LargePatches, col = "blue")
points(x = seq(-150, 150, by = 3), y = MedPatches, col = "darkred")
points(x = -150:150, y = SmallPatches, col = "black")

writeLines("Pass", con = "ModelTesting/EnvCalcFxns.txt")
