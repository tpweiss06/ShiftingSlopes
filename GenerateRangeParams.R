# This script will generate a few convenience functions to identify range parameters
#	for use in subsequent simulations. First values for delta, gamma, and local selection
#	will be arbitrarily set, and then values for tau will be calculated accordingly.
# The script will generate both a table with all the range parameter values to be used and
#	a visualization of the different parameter combinations.

setwd("~/Desktop/RangeShifts/ShiftingSlopesOther")

# Load the necessary libraries to put a table in an R graphic
library(grid)
library(gridBase)
library(gridExtra)

# First make a function to generate the realized patch carrying capacities throughout
#    the range based on the range parameters
K <- function(beta, gamma, tau, xVals, Kmax){
     f <- rep(NA, length(xVals))
     for(i in 1:length(xVals)){
          if(xVals[i] > beta){
               numerator <- exp(-1*gamma * (xVals[i] - beta - tau))
               denominator <- 1 + exp(-1*gamma * (xVals[i] - beta - tau))
               f[i] <- numerator / denominator
          } else if(xVals[i] <= beta){
               numerator <- exp(gamma * (xVals[i] - beta + tau))
               denominator <- 1 + exp(gamma * (xVals[i] - beta + tau))
               f[i] <- numerator / denominator
          } 
     }
     Kseq <- f * Kmax
     return(Kseq)
}

# Given a delta and minimum gamma value, calculate the tau value necessary for a
#    smooth function
CalcTauMin <- function(delta, gamma){
     numerator <- log( (1 - delta) / delta )
     TauMin <- numerator / gamma
     return(TauMin)
}

# Calculate a new tau value constrained by maintaining the width of the range
NewTau <- function(xHat, gamma, fHat){
     TauNew <- xHat + (1/gamma) * log(fHat/(1-fHat))
     return(TauNew)
}

# Calculate the total area under a given range curve
Integral <- function(gamma, tau){
     numerator <- 2 * log( exp(gamma * tau) + 1 )
     Area <- numerator / gamma
     return(Area)
}

# Finally, make a function to calculate a new Kmax to keep the total achievable
#    population size throughout the range constant.
NewK <- function(gamma, tau, TotalK){
     Area <- Integral(gamma = gamma, tau = tau)
     Kmax <- TotalK / Area
     return(Kmax)
}

# Now set the values for gamma, lambda, eta, and beta
GammaSeq <- c(0.0025, 0.0075, 0.25)
LambdaSeq <- c(0, 0.004, 0.008)
Beta <- 0
Eta <- 50

# Now create a sequence of x values and an object to store the tau values
xSeq <- seq(-2000, 2000, length.out = 10000) 
TauSeq <- rep(NA, length(GammaSeq))
# Set the initial tau value based on previous exploratory trials
TauSeq[1] <- -240

# Now calculate the values for ranges that keep the width constant. First, define
#    the width of the range as the interval in which f(x,t) is equal or greater than
#    0.1 (i.e. 10% of the maximum patch carrying capacity). Now use this definition
#    to find the x values that will correspond to this threshold in all ranges.
fHat <- 0.1
xHat <- TauSeq[1] - (1/GammaSeq[1])*log(fHat/(1-fHat))

# Now calculate the rest of the tau values corresponding to this width
for(i in 2:length(GammaSeq)){
     TauSeq[i] <- xHat + (1/GammaSeq[i])*log(fHat/(1-fHat))
}

# Using these tau values, calculate the area under the f(x,t) curves and use that
#    to calculate the different Kmax values for each scenario
AreaSeq <- rep(NA, length(GammaSeq))
for(i in 1:3){
     AreaSeq[i] <- Integral(gamma = GammaSeq[i], tau = TauSeq[i])
}
# Now arbitrarily set the Kmax for the gradual range edge and use this to set the
#    others. This value is based on previous exploratory trials
KmaxSeq <- rep(NA, length(GammaSeq))
KmaxSeq[1] <- 240
TotalK <- AreaSeq[1] * KmaxSeq[1]
for(i in 2:length(GammaSeq)){
     KmaxSeq[i] <- NewK(gamma = GammaSeq[i], tau = TauSeq[i], TotalK = TotalK)
}

# Use these values to calculate the realized patch carrying capacities
#    throughout the ranges
PatchK <- matrix(NA, nrow = length(GammaSeq), ncol = length(xSeq))
for(i in 1:length(GammaSeq)){
     PatchK[i,] <- K(beta = Beta, gamma = GammaSeq[i], tau = TauSeq[i], 
                     xVals = xSeq, Kmax = KmaxSeq[i])
}

# Now, create a table to store the parameter values in
RangeParams <- expand.grid(gamma = GammaSeq, lambda = LambdaSeq)
RangeParams <- cbind(RangeParams, tau = rep(TauSeq, 3), Kmax = rep(KmaxSeq, 3),
                     eta = rep(Eta, 9))
RangeParams
# Save these parameters for later use and put them in a table grob for plotting
write.csv(RangeParams, file = "RangeParameters.csv")
PlotTable <- RangeParams
PlotTable[,3] <- round(PlotTable[,3], digits = 3)
PlotTable[,4] <- round(PlotTable[,4], digits = 4)
PlotTableGrob <- tableGrob(d = PlotTable, cols = NULL)

# Finally, calculate the environmental niche values throughout the range
PatchBounds <- range(xSeq) / Eta
PatchCenter <- seq(PatchBounds[1], PatchBounds[2])
Zopt <- matrix(NA, nrow = length(GammaSeq), ncol = length(PatchCenter))
for(i in 1:length(GammaSeq)){
     Zopt[i,] <- LambdaSeq[i] * (PatchCenter * Eta - Beta)
}

# Now plot a visualization of the range parameters
pdf(file = "SchematicFigures/RangeParams.pdf", width = 7, height = 4, onefile = FALSE, paper = "special")
     par(mfrow = c(1,2), oma = c(0, 1, 0, 0), mar = c(5, 4.25, 4, 2.25) + 0.1)
     # Make the K figure
     plot(x = xSeq, y = PatchK[1,], type = "l", col = "red", lwd = 1.5, main = "", xlab = "",
          ylab = "K", xaxt = "n", las = 1, cex.axis = 1.5, ylim = c(0,120), cex.lab = 1.25)
     lines(x = xSeq, y = PatchK[2,], lwd = 1.5, col = "blue")
     lines(x = xSeq, y = PatchK[3,], lwd = 1.5, col = "green")
     abline(v = Beta, lty = 2, lwd = 0.75)
     text(x = 0, y = -15, labels = expression(beta["t"]), cex = 1.5, xpd = NA)
     
     # Make the Zopt figure
     plot(x = PatchCenter, y = Zopt[3,], type = "l", col = "black", main = "", xlab = "",
          ylab = expression("Z"["opt"]), xaxt = "n", las = 1, cex.axis = 1.5, cex.lab = 1.25)
     lines(x = PatchCenter, y = Zopt[1,])
     lines(x = PatchCenter, y = Zopt[2,])
     abline(v = Beta, lty = 2, lwd = 0.75)
     text(x = 0, y = -20, labels = expression(beta["t"]), cex = 1.5, xpd = NA)
     
     # Now add in the x label for both plots
     mtext("Spatial location (x)", side = 1, outer = TRUE, cex = 1.5, line = -2)
     #text(x = range(xSeq)[1] * 2, y = -10, labels = "Spatial location (x)", 
     #     xpd = NA, cex = 1.25)
     
     # Finally, add the range parameters
     # Use frame() to move on to the next graphing window
     #frame()
     #vps <- baseViewports()
     #pushViewport(vps$inner, vps$figure, vps$plot)
     #grid.draw(PlotTableGrob)
     #text(x = 0.15, y = 1.1, labels = expression(gamma), xpd = NA, cex = 2)
     #text(x = 0.6, y = 1.1, labels = expression(lambda), xpd = NA, cex = 2)
     #text(x = 1.05, y = 1.1, labels = expression(tau), xpd = NA, cex = 2)
     #text(x = 1.5, y = 1.1, labels = expression("K"["max"]))
     #text(x = 2, y = 1.1, labels = expression(eta), xpd = NA, cex = 2)
dev.off()

# Now create a second schematic figure to illustrate the scale of discretization
#    in the landscapes
pdf(file = "SchematicFigures/Discretization.pdf", width = 7, height = 7, onefile = FALSE, paper = "special")
     par(mfrow = c(3,1), mar = c(3,5,3,2) + 0.1, oma = c(2,0,0,0))
     # Make the f(x) figures
     for(i in 1:3){
          #CenterVals <- f_xt(beta = Beta, gamma = GammaSeq[i], tau = TauSeq[i], xVals = PatchCenter*Eta)
          CenterVals <- K(beta = Beta, gamma = GammaSeq[i], tau = TauSeq[i], 
                          xVals = PatchCenter*Eta, Kmax = KmaxSeq[i])
          
          plot(NA, NA, xlim = range(xSeq), ylim = c(0, 120), main = "", xlab = "",
               ylab = "K", las = 1, cex.lab = 1.25) #xaxt = "n"
          points(x = PatchCenter*Eta, y = CenterVals, col = "black", pch = 16)
          abline(v = xHat, lty = 2)
          abline(v = -1*xHat, lty = 2)
          
          if(i == 1){
               mtext(expression(paste(gamma, " = 0.0025", sep = "")), side = 3,
                     line = 0, cex = 1.25)
          } else if(i == 2){
               mtext(expression(paste(gamma, " = 0.0075", sep = "")), side = 3,
                     line = 0, cex = 1.25)
          } else{
               mtext(expression(paste(gamma, " = 0.25", sep = "")), side = 3,
                     line = 0, cex = 1.25)
               mtext("Spatial location (x)", side = 1, cex = 1.5, line = 3)
          }
     }
dev.off()

