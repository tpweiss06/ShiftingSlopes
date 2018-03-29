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

# Make a function to calculate the range curve over a given spatial extent
f_xt <- function(beta, gamma, tau, xVals){
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
     return(f)
}

# Given a delta and minimum gamma value, calculate the tau value necessary for a
#    smooth function
CalcTauMin <- function(delta, gamma){
     numerator <- log( (1 - delta) / delta )
     TauMin <- numerator / gamma
     return(TauMin)
}

# When changing the slope of the range boundary, keep the area of the range
#    constant and calculate a new tau
CalcTauNew <- function(gamma, C){
     numerator <- log( exp( (C*gamma) / 2 ) - 1 )
     TauNew <- numerator / gamma
     return(TauNew)
}

# Calculate the total area under a given range curve
Calc_C <- function(gamma, tau){
     numerator <- 2 * log( exp(gamma * tau) + 1 )
     C <- numerator / gamma
     return(C)
}


# Now set the values for gamma, delta, eta, and beta
GammaSeq <- c(0.0025, 0.025, 0.25)
Delta <- 1e-3
Beta <- 0
Eta <- 50

# Now create a sequence of x values and an object to store the tau values
xSeq <- seq(-2000, 2000, length.out = 10000) 
TauSeq <- rep(NA, length(GammaSeq))

# Instead of calculating an initial tau according to Delta, I'm going to arbitrarily
#    pick a value that won't produce a smooth curve and go from there. This will
#    allow an exploration of a wider range of range edge conditions
TauSeq[1] <- 250
C <- Calc_C(gamma = GammaSeq[1], tau = TauSeq[1])
TauSeq[2] <- CalcTauNew(gamma = GammaSeq[2], C = C)
TauSeq[3] <- CalcTauNew(gamma = GammaSeq[3], C = C)

# Now create a vector for the lambda values to use for local selection
#LambdaSeq <- c(0, 0.01, 0.02)
LambdaSeq <- c(0, 0.004, 0.008)

# Now, create a table to store the parameter values in
RangeParams <- expand.grid(gamma = GammaSeq, lambda = LambdaSeq)
RangeParams <- cbind(RangeParams, tau = rep(TauSeq, 3), eta = rep(Eta, 9))
RangeParams
# Save these parameters for later use and put them in a table grob for plotting
write.csv(RangeParams, file = "RangeParameters.csv")
PlotTable <- RangeParams
PlotTable[,3] <- round(PlotTable[,3], digits = 3)
PlotTableGrob <- tableGrob(d = PlotTable, cols = NULL)

# Finally, use these values to calculate starting range shapes for each
#    slope value
PatchBounds <- range(xSeq) / Eta
PatchCenter <- seq(PatchBounds[1], PatchBounds[2])

fseq <- matrix(NA, nrow = length(GammaSeq), ncol = length(xSeq))
Zopt <- matrix(NA, nrow = length(GammaSeq), ncol = length(PatchCenter))
for(i in 1:length(GammaSeq)){
     fseq[i,] <- f_xt(beta = Beta, gamma = GammaSeq[i], tau = TauSeq[i], xVals = xSeq)
     Zopt[i,] <- LambdaSeq[i] * (PatchCenter * Eta - Beta)
}

# Now plot a visualization of the parameters used
# First make a layout matrix for the graph
FigMat <- matrix(NA, nrow = 2, ncol = 6)
FigMat[1,] <- c(1,1,2,2,3,3)
FigMat[2,] <- c(4,4,4,5,5,6)

pdf(file = "SchematicFigures/RangeParams.pdf", width = 7, height = 7, onefile = FALSE, paper = "special")
     par(oma = c(0,2,0,0))
     layout(mat = FigMat)
     # Make the f(x) figures
     for(i in 1:3){
          plot(x = xSeq, y = fseq[i,], type = "l", col = "black", main = "", xlab = "",
               ylab = "", xaxt = "n", las = 1, cex.axis = 1.5, ylim = c(0,1))
          # Add vertical lines in for beta and beta +- tau with labels
          segments(x0 = Beta, y0 = -1, x1 = Beta, y1 = 1, lty = 2)
          text(x = 0, y = -0.1, labels = expression(beta["t"]), cex = 1.5, xpd = NA)
          
          if(i == 1){
               mtext("f(x,t)", side = 2, cex = 1.5, line = 3.5)
               mtext(expression(paste(gamma, " = 0.0025", sep = "")), side = 3,
                     line = 0, cex = 1.25)
          } else if(i == 2){
               mtext("Spatial location (x)", side = 1, cex = 1.5, line = 4)
               mtext(expression(paste(gamma, " = 0.025", sep = "")), side = 3,
                     line = 0, cex = 1.25)
          } else{
               mtext(expression(paste(gamma, " = 0.25", sep = "")), side = 3,
                     line = 0, cex = 1.25)
          }
     }
     
     # Make the Zopt figure
     PlotBottom <- -50
     PlotTop <- 50
     plot(x = PatchCenter, y = Zopt[3,], type = "l", col = "black", main = "", xlab = "",
          ylab = "", xaxt = "n", las = 1, cex.axis = 1.5)
     mtext(expression("Z"["opt"]), side =2, line = 2.5, cex = 1.5)
     lines(x = PatchCenter, y = Zopt[1,])
     lines(x = PatchCenter, y = Zopt[2,])
     # Add vertical lines in for beta
     segments(x0 = Beta, y0 = PlotBottom + 5, x1 = Beta, y1 = PlotTop, lty = 2)
     text(x = 0, y = PlotBottom - 5, labels = expression(beta["t"]), cex = 1.5)
     
     # Add in the lambda values
     text(x = 3800, y = 45, labels = expression(paste(lambda, " = 0.01", sep = "")), 
          cex = 1.25, srt = 40)
     text(x = 4200, y = 27, labels = expression(paste(lambda, " = 0.005", sep = "")), 
          cex = 1.25, srt = 25)
     text(x = 3800, y = 5, labels = expression(paste(lambda, " = 0", sep = "")), 
          cex = 1.25)
     mtext("Spatial location", side = 1, line = 1, cex = 1)
     
     # Finally, add the range parameters
     # Use frame() to move on to the next graphing window
     frame()
     vps <- baseViewports()
     pushViewport(vps$inner, vps$figure, vps$plot)
     grid.draw(PlotTableGrob)
     text(x = 0.15, y = 1.1, labels = expression(gamma), xpd = NA, cex = 2)
     text(x = 0.6, y = 1.1, labels = expression(lambda), xpd = NA, cex = 2)
     text(x = 1.05, y = 1.1, labels = expression(tau), xpd = NA, cex = 2)
     text(x = 1.5, y = 1.1, labels = expression(eta), xpd = NA, cex = 2)
dev.off()

# Now create a second schematic figure to illustrate the scale of discretization
#    in the landscapes
pdf(file = "SchematicFigures/Discretization.pdf", width = 7, height = 7, onefile = FALSE, paper = "special")
     par(mfrow = c(3,1), mar = c(3,5,3,2) + 0.1, oma = c(2,0,0,0))
     # Make the f(x) figures
     for(i in 1:3){
          CenterVals <- f_xt(beta = Beta, gamma = GammaSeq[i], tau = TauSeq[i], xVals = PatchCenter*Eta)
          ColVals <- ifelse(CenterVals < 0.02, "black", "red")
          plot(NA, NA, xlim = range(xSeq), ylim = c(0, 1), main = "", xlab = "",
               ylab = "f(x,t)", las = 1) #xaxt = "n"
          points(x = PatchCenter*Eta, y = CenterVals, col = ColVals, pch = 16)
          
          legend("topright", bty = "n", legend = c("Viable Patch", "Unviable Patch"),
                 col = c("red", "black"), pch = 16)
          
          if(i == 1){
               mtext(expression(paste(gamma, " = 0.0025", sep = "")), side = 3,
                     line = 0, cex = 1.25)
          } else if(i == 2){
               mtext(expression(paste(gamma, " = 0.025", sep = "")), side = 3,
                     line = 0, cex = 1.25)
          } else{
               mtext(expression(paste(gamma, " = 0.25", sep = "")), side = 3,
                     line = 0, cex = 1.25)
               mtext("Spatial location (x)", side = 1, cex = 1.5, line = 3)
          }
     }
dev.off()

