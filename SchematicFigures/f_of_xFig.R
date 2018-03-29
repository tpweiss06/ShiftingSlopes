# This script makes a nice, annnotated graph of f(x) with the relevant parameters
#    labeled.

setwd("~/Desktop/RangeShifts/ShiftingSlopes/SchematicFigures/")

# Make a function to calculate f(x) over a given spatial extent
RangeCapacity <- function(beta, gamma, tau, xSeq){
     f <- rep(NA, length(xSeq))
     for(i in 1:length(xSeq)){
          if(xSeq[i] > beta){
               numerator <- exp(-1*gamma * (xSeq[i] - beta - tau))
               denominator <- 1 + exp(-1*gamma * (xSeq[i] - beta - tau))
               f[i] <- numerator / denominator
          } else if(xSeq[i] <= beta){
               numerator <- exp(gamma * (xSeq[i] - beta + tau))
               denominator <- 1 + exp(gamma * (xSeq[i] - beta + tau))
               f[i] <- numerator / denominator
          } 
     }
     return(f)
}

# Now pick some parameters to use for the visualization and calculate f(x)
xSeq <- seq(-500, 500, length.out = 10000)
Beta <- 0
Gamma <- 0.03
Tau <- 250
Eta <- 25
f <- RangeCapacity(beta = Beta, gamma = Gamma, tau = Tau, xSeq = xSeq)
FigMat <- matrix(c(1,1,2), nrow = 3, ncol = 1)

# Now make a nice plot of f(x) with parameter values annotated
pdf(file = "f_of_xt.pdf", width = 7, height = 3.5, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = c(0.5, 7, 0.5, 2) + 0.1, oma = c(3,0,0,0))
     plot(x = xSeq, y = f, type = "l", col = "black", main = "", xlab = "",
          ylab = "", xaxt = "n", las = 1, cex.axis = 1.5)
     # Add a vertical line in for beta with a label
     segments(x0 = Beta, y0 = -1, x1 = Beta, y1 = 1, lty = 2)
     #abline(v = Beta, lty = 2)
     text(x = 20, y = 0.02, labels = expression(beta["t"]), cex = 1.5)
     # Add horizontal lines in for beta +- tau
     segments(x0 = Beta + Tau, y0 = -1, x1 = Beta + Tau, y1 = 0.5, lty = 2)
     segments(x0 = Beta - Tau, y0 = -1, x1 = Beta - Tau, y1 = 0.5, lty = 2)
     text(x = Beta+Tau - 40, y = 0.02, labels = expression(beta["t"] + tau), cex = 1.5,
          xpd = NA)
     text(x = Beta-Tau + 40, y = 0.02, labels = expression(beta["t"] - tau), cex = 1.5,
          xpd = NA)
     # Add tangent lines for the slope at beta +- tau (or possibly f'(x) labels)
     # First the Beta + Tau side
     delta <- 60
     m <- (-1 * Gamma) / 4
     xHat <- Beta + Tau
     yHat <- 0.5
     y1 <- yHat + delta*m
     y2 <- yHat - delta*m
     segments(x0 = xHat + delta, y0 = y1, x1 = xHat - delta, y1 = y2, lty = 2)
     text(x = xHat + delta, y = 0.5, labels = expression(-gamma / 4), cex = 1.5)
     
     # Now the Beta - Tau side
     xHat <- Beta - Tau
     m <- Gamma / 4
     y1 <- yHat + delta*m
     y2 <- yHat - delta*m
     segments(x0 = xHat + delta, y0 = y1, x1 = xHat - delta, y1 = y2, lty = 2)
     text(x = xHat - delta, y = 0.5, labels = expression(gamma / 4), cex = 1.5)
     
     # Axis label
     mtext(expression(italic("f(x,t)")), side = 2, line = 3, cex = 1.5, las = 2)
     
     
     # Now add the grid for the discretized landscape
     plot(NA, NA, xlim = range(xSeq), ylim = c(-1*Eta, 6*Eta), axes = FALSE,
          xlab = "", ylab = "")
     # Put in the horizontal segments
     for(i in 0:5){
          segments(x0 = -400 - Eta/2, y0 = i*Eta - Eta/2, x1 = 400 + Eta/2, y1 = i*Eta - Eta/2)
     }
     # Now put in the vertical segments
     LeftEdge <- seq(-412.5, 387.5, by = Eta)
     RightEdge <- seq(-387.5, 412.5, by = Eta)
     for(i in 1:length(LeftEdge)){
          segments(x0 = LeftEdge[i], y0 = 0 - Eta/2, x1 = LeftEdge[i], y1 = 5*Eta - Eta/2)
          segments(x0 = RightEdge[i], y0 = 0 - Eta/2, x1 = RightEdge[i], y1 = 5*Eta - Eta/2)
     }
     text(labels = "Discretized\nLandscape", x = -600, y = 50, xpd = NA, cex = 2)
     
     # Create an enlarged patch for visualization
     polygon(x = c(400 + 2*Eta, 400 + 2*Eta, 400 + 5*Eta, 400 + 5*Eta),
             y = c(2*Eta, 5*Eta, 5*Eta, 2*Eta), xpd = NA)
     segments(x0 = 400 - Eta/2, y0 = 3*Eta - Eta/2, x1 = 400 + 2*Eta, y1 = 5*Eta)
     segments(x0 = 400 - Eta/2, y0 = 2*Eta - Eta/2, x1 = 400 + 2*Eta, y1 = 2*Eta)
     segments(x0 = 400 + Eta/2, y0 = 2*Eta - Eta/2, x1 = 400 + 5*Eta, y1 = 2*Eta)
     text(labels = expression(eta), x = 400+3.5*Eta, y = 5*Eta+20, xpd = NA, cex = 1.5)
     text(labels = expression(eta), x = 400+5*Eta+15, y = 3.5*Eta, xpd = NA, cex = 1.5)
     
     # Now add arrows to illustrate dispersal and wrapping boundaries
     arrows(x0 = -400, y0 = 0, x1 = -400 - 1.25*Eta, y1 = 0, length = 0.05)
     arrows(x0 = 400, y0 = 0, x1 = 400 + 1.25*Eta, y1 = 0, length = 0.05)
     arrows(x0 = 0, y0 = -1.05*Eta, x1 = 0, y1 = 0 + 0.05*Eta, length = 0.05, code = 3)
     arrows(x0 = 0, y0 = 3.95*Eta, x1 = 0, y1 = 5.05*Eta, length = 0.05, code = 3)
     
     mtext("Space (x)", side = 1, line = 2, cex = 1.5)
dev.off()

