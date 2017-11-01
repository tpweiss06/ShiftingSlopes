# This script will generate a few convenience functions to identify range parameters
#	for use in subsequent simulations. First values for delta, gamma, and local selection
#	will be arbitrarily set, and then values for tau will be calculated accordingly.
# The script will generate both a table with all the range parameter values to be used and
#	a visualization of the different parameter combinations.
# Note: For simplicity, alpha is assumed to be 1 here and is therefore removed from the following
#	functions.

setwd("~/Desktop/RangeShifts/ShiftingSlopes")

# Make a function to calculate the range curve over a given spatial extent
RangeCapacity <- function(beta, gamma, tau, x_vals){
     R_seq <- rep(NA, length(x_vals))
     for(i in 1:length(x_vals)){
          if(x_vals[i] > beta){
               numerator <- exp(-1*gamma * (x_vals[i] - beta - tau))
               denominator <- 1 + exp(-1*gamma * (x_vals[i] - beta - tau))
               R_seq[i] <- numerator / denominator
          } else if(x_vals[i] < beta){
               numerator <- exp(gamma * (x_vals[i] - beta + tau))
               denominator <- 1 + exp(gamma * (x_vals[i] - beta + tau))
               R_seq[i] <- numerator / denominator
          } else if(x_vals[i] == beta){
               R_seq[i] <- 1
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
Cubic <- function(beta, tau, x_vals, LocalSel){
     Zopt <- LocalSel * (x_vals / (0.5*tau)) ^ 3
     #Zopt <- LocalSel * (x_vals / (tau)) ^ 3
     return(Zopt)
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


# Now set the values for gamma, delta, beta, and LocalSel
GammaSeq <- c(0.0025, 0.01, 0.25)
LocalSel <- c(0.25, 1, 1.75)
Delta <- 1e-3
Beta <- 0

# Now create a sequence of x values and an object to store the tau values
x_seq <- seq(-5000, 5000, length.out = 10000) 
TauSeq <- rep(NA, length(GammaSeq))

# Now calculate TauMin and C
TauSeq[1] <- CalcTauMin(delta = Delta, gamma = GammaSeq[1])
C <- Calc_C(gamma = GammaSeq[1], tau = TauSeq[1])
TauSeq[2] <- CalcTauNew(gamma = GammaSeq[2], C = C)
TauSeq[3] <- CalcTauNew(gamma = GammaSeq[3], C = C)

# Now, create a table to store the parameter values in
RangeParams <- expand.grid(gamma = GammaSeq, LocalSel = LocalSel)
RangeParams <- cbind(RangeParams, tau = rep(TauSeq, 3))
RangeParams
# Save these parameters for later use
write.csv(RangeParams, file = "RangeParameters.csv")


# Finally, use these values to calculate starting range shapes for each
#    slope value
Rseq <- matrix(NA, nrow = length(GammaSeq), ncol = length(x_seq))
FitVals <- matrix(NA, nrow = length(GammaSeq), ncol = length(x_seq))
CubicVals <- matrix(NA, nrow = length(GammaSeq), ncol = length(x_seq))
for(i in 1:length(GammaSeq)){
     Rseq[i,] <- RangeCapacity(beta = Beta, gamma = GammaSeq[i], 
                        tau = TauSeq[i], x_vals = x_seq)
     FitVals[i,] <- RelFit(beta = Beta, tau = TauSeq[i], x_vals = x_seq, LocalSel = LocalSel[i])
     CubicVals[i,] <- Cubic(beta = Beta, tau = TauSeq[i], x_vals = x_seq, LocalSel = LocalSel[i])
}

# Now plot a visualization of the parameters used
pdf(file = "RangeParams.pdf", width = 7, height = 7, onefile = FALSE, paper = "special")
	par(mfrow = c(3,3), mar = c(0,0,0,0), oma = c(4, 4, 4,6))
	for(i in 1:3){
		for(j in 1:3){
			plot(x_seq, Rseq[j,], type = "l", axes = FALSE, col = "blue")
			lines(x_seq, FitVals[i,], col = "red")
			box()
			if(j == 3){
				LambdaText <- bquote(lambda == .(LocalSel[i]))
				mtext(LambdaText, side = 4, col = "red", cex = 1.25, line = 1.5)
			}
			if(i == 1){
				GammaText <- bquote(gamma == .(GammaSeq[j]))
				mtext(GammaText, side = 3, col = "blue", cex = 1.25)
			}
		}
	}
	mtext("Spatial Position", outer = TRUE, side = 1, cex = 1.5, line = 1.5)
	mtext("Range Capacity", outer = TRUE, side = 2, cex = 1.5, col = "blue", line = 1.5)
	mtext("Local Optimum", outer = TRUE, side = 4, cex = 1.5, col = "red", line = 3.5)
dev.off()

# Now plot a visualization of the parameters used with the cubic Zopt function
pdf(file = "RangeParamsCubic1.pdf", width = 7, height = 7, onefile = FALSE, paper = "special")
     par(mfrow = c(3,3), mar = c(0,0,0,0), oma = c(4, 4, 4,6))
     for(i in 1:3){
          for(j in 1:3){
               plot(x_seq, Rseq[j,], type = "l", axes = FALSE, col = "blue")
               par(new = TRUE)
               plot(x_seq, CubicVals[i,], col = "red", axes = FALSE, type = "l",
                    xlim = range(x_seq), ylim = c(-11,11))
               box()
               if(j == 3){
                    LambdaText <- bquote(lambda == .(LocalSel[i]))
                    mtext(LambdaText, side = 4, col = "red", cex = 1.25, line = 1.5)
               }
               if(i == 1){
                    GammaText <- bquote(gamma == .(GammaSeq[j]))
                    mtext(GammaText, side = 3, col = "blue", cex = 1.25)
               }
          }
     }
     mtext("Spatial Position", outer = TRUE, side = 1, cex = 1.5, line = 1.5)
     mtext("Range Capacity", outer = TRUE, side = 2, cex = 1.5, col = "blue", line = 1.5)
     mtext("Local Optimum", outer = TRUE, side = 4, cex = 1.5, col = "red", line = 3.5)
dev.off()


