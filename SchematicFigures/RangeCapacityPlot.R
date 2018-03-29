setwd("~/Desktop/RangeShifts/")

# Make a function to calculate the range curve over a given spatial extent
RangeCapacity <- function(alpha, beta, gamma, tau, x_vals){
     R_seq <- rep(NA, length(x_vals))
     for(i in 1:length(x_vals)){
          if(x_vals[i] > beta){
               numerator <- alpha * exp(-1*gamma * (x_vals[i] - beta - tau))
               denominator <- 1 + exp(-1*gamma * (x_vals[i] - beta - tau))
               R_seq[i] <- numerator / denominator
          } else if(x_vals[i] < beta){
               numerator <- alpha * exp(gamma * (x_vals[i] - beta + tau))
               denominator <- 1 + exp(gamma * (x_vals[i] - beta + tau))
               R_seq[i] <- numerator / denominator
          } else if(x_vals[i] == beta){
               R_seq[i] <- alpha
          }
     }
     return(R_seq)
}

# Make a function to display the local fitness optimum across the range (lower
#    values don't mean lower fitness, just a different optimum phenotype; actual
#    fitness values will be calculated with these local optima and stabilizing
#    selection; reproduction will occur according to the range capacity gradient)
RelFit <- function(beta, tau, x_vals, LocalSel){
     FitSeq <- rep(NA, length(x_vals))
     GammaVal <- LocalSel / tau
     for(i in 1:length(x_vals)){
          numerator <- exp(GammaVal * (x_vals[i] - beta))
          denominator <- 1 + exp(GammaVal * (x_vals[i] - beta))
          FitSeq[i] <- numerator / denominator
     }
     return(FitSeq)
}


# Given a delta and minimum gamma value, calculate the tau value necessary for a
#    smooth function
CalcTauMin <- function(alpha, delta, gamma){
     numerator <- log( (alpha - delta) / delta )
     TauMin <- numerator / gamma
     return(TauMin)
}

# When changing the slope of the range boundary, keep the area of the range
#    constant and calculate a new tau
CalcTauNew <- function(alpha, gamma, C){
     numerator <- log( exp( (C*gamma) / (2*alpha) ) - 1 )
     TauNew <- numerator / gamma
     return(TauNew)
}

# Calculate the total area under a given range curve
Calc_C <- function(alpha, gamma, tau){
     numerator <- 2 * alpha * log( exp(gamma * tau) + 1 )
     C <- numerator / gamma
     return(C)
}


# Now pick some parameters to use for the initial visualization
x_seq <- seq(-500, 500, length.out = 10000)
Alpha <- 1
Delta <- 1e-3
BetaSeq <- seq(0, 500, by = 100)
GammaSeq <- c(0.025, 0.075, 0.25)
TauSeq <- rep(NA, length(GammaSeq))

# Now calculate TauMin and C
TauSeq[1] <- CalcTauMin(alpha = Alpha, delta = Delta, gamma = GammaSeq[1])
C <- Calc_C(alpha = Alpha, gamma = GammaSeq[1], tau = TauSeq[1])

# Now calculate the other tau values
for(i in 2:length(TauSeq)){
     TauSeq[i] <- CalcTauNew(alpha = Alpha, C = C, gamma = GammaSeq[i])
}

# Finally, use these values to calculate a starting range shapes for each
#    slope value
Rseq <- matrix(NA, nrow = length(GammaSeq), ncol = length(x_seq))
FitVals <- matrix(NA, nrow = length(GammaSeq), ncol = length(x_seq))
for(i in 1:length(GammaSeq)){
     Rseq[i,] <- RangeCapacity(alpha = Alpha, beta = BetaSeq[1], gamma = GammaSeq[i], 
                        tau = TauSeq[i], x_vals = x_seq)
     FitVals[i,] <- RelFit(beta = BetaSeq[1], tau = TauSeq[i], x_vals = x_seq, LocalSel = 0.25)
}

# Try a basic plot to make sure these functions are working and that these
#    initial values work
plot(x_seq, Rseq[1,], type = "l")
lines(x_seq, FitVals[1,])
plot(x_seq, Rseq[2,], type = "l")
lines(x_seq, FitVals[2,])
plot(x_seq, Rseq[3,], type = "l")
lines(x_seq, FitVals[3,])

# Now make a nice plot of the stationary range capacity and a visualization of 
#    the shift itself
PlotMatrix <- matrix(c(1,2,2,3,4,4,5,6,6), nrow = 3, ncol = 3, byrow = TRUE)
#layout(mat = PlotMatrix)
#layout.show(n=6)
TimeCols <- rainbow(n = length(BetaSeq))
newX <- seq(-500, 1000, length.out = 1000000)
RangeCol <- "lightcoral"

pdf(file = "RangeCapacity.pdf", width = 7, height = 5, onefile = FALSE, paper = "special")
     layout(mat = PlotMatrix)
     par(mar = c(1, 0.8, 0.8, 0.4) + 0.02, oma = c(5, 4, 4, 2))
     
     # First plot the shallowest slope in a stationary plot, then the shifting
     plot(x = x_seq, y = Rseq[1,], type = "l", lwd = 2, axes = FALSE, xlab = "",
          ylab = "", ylim = c(0, 1.1*Alpha))
     box()
     polygon(x = c(x_seq, rev(x_seq)), y = c(Rseq[1,], rep(0, length(x_seq))), 
             col = RangeCol, border = NA)
     mtext("Range Capacity", side = 3, cex = 1.5, line = 1)
     
     plot(x = NA, y = NA, axes = FALSE, xlab = "", ylab = "", xlim = range(newX),
          ylim = c(0, 1.1*Alpha))
     for(t in 1:length(BetaSeq)){
          Rnew <- RangeCapacity(alpha = Alpha, beta = BetaSeq[t], gamma = GammaSeq[1], 
                                tau = TauSeq[1], x_vals = newX)
          lines(x = newX, y = Rnew, col = TimeCols[t])
     }
     text(expression(italic("t = 0")), x = -400, y = 0.9*Alpha, cex = 1.25, col = TimeCols[1])
     text(expression(italic("t = 5")), x = 900, y = 0.9*Alpha, cex = 1.25, col = TimeCols[6])
     mtext("Range Shifting", side = 3, cex = 1.5, line = 1)
     box()

     # Now plot the second level
     plot(x = x_seq, y = Rseq[2,], type = "l", lwd = 2, axes = FALSE, xlab = "",
          ylab = "", ylim = c(0, 1.1*Alpha))
     box()
     polygon(x = c(x_seq, rev(x_seq)), y = c(Rseq[2,], rep(0, length(x_seq))), 
             col = RangeCol, border = NA)
     
     plot(x = NA, y = NA, axes = FALSE, xlab = "", ylab = "", xlim = range(newX),
          ylim = c(0, 1.1*Alpha))
     for(t in 1:length(BetaSeq)){
          Rnew <- RangeCapacity(alpha = Alpha, beta = BetaSeq[t], gamma = GammaSeq[2], 
                                tau = TauSeq[2], x_vals = newX)
          lines(x = newX, y = Rnew, col = TimeCols[t])
     }
     text(expression(italic("t = 0")), x = -400, y = 0.9*Alpha, cex = 1.25, col = TimeCols[1])
     text(expression(italic("t = 5")), x = 900, y = 0.9*Alpha, cex = 1.25, col = TimeCols[6])
     box()

     # Now plot the third level
     plot(x = x_seq, y = Rseq[3,], type = "l", lwd = 2, axes = FALSE, xlab = "",
          ylab = "", ylim = c(0, 1.1*Alpha))
     box()
     polygon(x = c(x_seq, rev(x_seq)), y = c(Rseq[3,], rep(0, length(x_seq))), 
             col = RangeCol, border = NA)
     
     plot(x = NA, y = NA, axes = FALSE, xlab = "", ylab = "", xlim = range(newX),
          ylim = c(0, 1.1*Alpha))
     for(t in 1:length(BetaSeq)){
          Rnew <- RangeCapacity(alpha = Alpha, beta = BetaSeq[t], gamma = GammaSeq[3], 
                                tau = TauSeq[3], x_vals = newX)
          lines(x = newX, y = Rnew, col = TimeCols[t])
     }
     text(expression(italic("t = 0")), x = -400, y = 0.9*Alpha, cex = 1.25, col = TimeCols[1])
     text(expression(italic("t = 5")), x = 900, y = 0.9*Alpha, cex = 1.25, col = TimeCols[6])
     box()
     
     mtext("Spatial location", side = 1, outer = TRUE, cex = 1.5, line = 2)
     mtext("Fitness metric", side = 2, outer = TRUE, cex = 1.5, line = 1)
dev.off()



