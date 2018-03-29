# This script will create the graphs from the abundance data.

# Set the working directory
setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/StationaryRange/")

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

# Set some graphical parameters to use for the subsequent figures
FigMat <- matrix(NA, nrow = 3, ncol = 13)
FigMat[1,] <- c(rep(1,4), rep(2,4), rep(3,4), 10)
FigMat[2,] <- c(rep(4,4), rep(5,4), rep(6,4), 10)
FigMat[3,] <- c(rep(7,4), rep(8,4), rep(9,4), 10)
OuterMar <- c(4, 10, 5, 2)
InnerMar <- c(1.5, 1.5, 1.5, 1.5)
TextSize <- 1.5
AxisSize <- 1.05
ArrowWidth <- 3
ArrowLength <- 0.25
FigWidth <- 8
FigHeight <- 6
LocLabels <- seq(-75, 75, by = 30)
maxt <- 2000
SimSeq <- c(7,8,9,4,5,6,1,2,3)
ColKeyWidth <- 15
xLabLine <- 1.75
yLabLine <- 1.5
GradLabLine <- 3
GradSubLine <- 1.5
AdaptLabLine <- 7
AdaptSubLine <- 5
TopArrow <- matrix(c(-3000, 1.2, 19000, 1.2), nrow = 2, ncol = 2, byrow = TRUE)
SideArrow <- matrix(c(-6200, -2.5, -6200, 0.8), nrow = 2, ncol = 2, byrow = TRUE)
LowAdj <- 0.05
HighAdj <- 0.95

# Set the patch sequence to use and use it later to set the continuous scale for
#    f(x,t)
PatchSeq <- seq(-60, 60, by = 1)

# Also get beta, tau, etc. from the range parameter csv file.
Beta <- 0
RangeParams <- read.csv("~/Desktop/RangeShifts/ShiftingSlopesOther/RangeParameters.csv")
xfseq <- seq(PatchSeq[1]*RangeParams$eta[1], 
             PatchSeq[length(PatchSeq)]*RangeParams$eta[1], length.out = 10000)

# Get the maximum values for the abundance values 
load("StationaryAbundResults.rdata")
MaxAbunds <- rep(NA, 9)
for(i in 1:9){
     MaxAbunds[i] <- max(SectorMean[i,,], na.rm = TRUE)
}

# Use the maxes determnined below to plot f(x,t) against patch_abund/max on a 0 to 1 scale
pdf(file = "Figures/Abund_fxt.pdf", width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          # plot the f(x,t) line
          fSeq <- RangeCapacity(beta = Beta, gamma = RangeParams$gamma[i], 
                                tau = RangeParams$tau[i], xSeq = xfseq)
          plot(xfseq, fSeq, xlim = range(xfseq), ylim = c(0,1), ylab = "", xlab = "",
               main = "", xaxt = "n", cex.axis = AxisSize, las = 1, type = "l",
               col = "blue")
          
          # add the points for realized abundances (scaled relative to the maximum
          #    realized abundance)
          AbundSeq <- SectorMean[i,,maxt]
          ScaledAbunds <- AbundSeq / MaxAbunds[i]
          PatchPoints <- seq(PatchSeq[1]*RangeParams$eta[1],
                             PatchSeq[length(PatchSeq)]*RangeParams$eta[1],
                             by = RangeParams$eta[1])
          points(x = PatchPoints, y = ScaledAbunds, pch = 19, col = "red")
          
          # Add the x-axis
          AxisTicks <- seq(PatchSeq[1]*RangeParams$eta[1],
                             PatchSeq[length(PatchSeq)]*RangeParams$eta[1],
                             by = RangeParams$eta[1]*10)
          Axislabels <- seq(PatchSeq[1], PatchSeq[length(PatchSeq)],
                            by = 10)
          axis(1, at = AxisTicks, labels = Axislabels, cex.axis = AxisSize)
          
          # Add the selection and environmental gradient arrows
          if(i == 7){
               arrows(x0 = TopArrow[1,1], y0 = TopArrow[1,2], x1 = TopArrow[2,1], 
                      y1 = TopArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
               arrows(x0 = SideArrow[1,1], y0 = SideArrow[1,2], x1 = SideArrow[2,1], 
                      y1 = SideArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
          }
     }

     # Add the x and y axis labels
     mtext("Spatial location (x)", side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
     mtext("Scaled Abundance", side = 2, outer = TRUE, line = yLabLine, cex = TextSize)

     # Add the selection and environmental gradient text
     mtext("Strength of local adaptation", side = 2, outer = TRUE, line = AdaptLabLine,
           cex = TextSize)
     mtext("Weak", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize,
           adj = LowAdj)
     mtext("Strong", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize,
           adj = HighAdj)
     
     mtext("Gradient at range edge", side = 3, outer = TRUE, line = GradLabLine,
           cex = TextSize)
     mtext("Gradual", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
           adj = LowAdj)
     mtext("Severe", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
           adj = HighAdj)
dev.off()


