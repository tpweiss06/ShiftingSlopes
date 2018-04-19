# This script will create the graphs from the extinction data.
library(RColorBrewer)

# Set the working directory
setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/ShiftingRange/")
source("~/Desktop/RangeShifts/ShiftingSlopesCode/ShiftingRange/ShiftParams.R")
RangeParams <- read.csv("~/Desktop/RangeShifts/ShiftingSlopesOther/RangeParameters.csv")

# Load the extinction data
load("Extinctions.rdata")
slow <- Extinctions[1,,]
med <- Extinctions[2,,]
fast <- Extinctions[3,,]
for(i in 1:9){
     slow[i,] <- cumsum(slow[i,]) / 100
     med[i,] <- cumsum(med[i,]) / 100
     fast[i,] <- cumsum(fast[i,]) / 100
}

# Set some graphical parameters to use for the subsequent figures
FigMat <- matrix(NA, nrow = 3, ncol = 12)
FigMat[1,] <- c(rep(1,4), rep(2,4), rep(3,4))
FigMat[2,] <- c(rep(4,4), rep(5,4), rep(6,4))
FigMat[3,] <- c(rep(7,4), rep(8,4), rep(9,4))
OuterMar <- c(4, 10, 5, 2)
InnerMar <- c(1.5, 1.5, 1.5, 1.5)
TextSize <- 1.5
AxisSize <- 1.05
ArrowWidth <- 3
ArrowLength <- 0.25
FigWidth <- 8
FigHeight <- 6
LocLabels <- seq(-60, 60, by = 24)
TimeSeq <- 1:100
TimeLabels <- seq(0, 100, by = 20)
SimSeq <- c(7,8,9,4,5,6,1,2,3)
ColKeyWidth <- 15
xLabLine <- 1.75
yLabLine <- 1.5
GradLabLine <- 3
GradSubLine <- 1.5
AdaptLabLine <- 7
AdaptSubLine <- 5
TopArrow <- matrix(c(10, 1.2, 740, 1.2), nrow = 2, ncol = 2, byrow = TRUE)
SideArrow <- matrix(c(-95, -2.75, -95, 0.9), nrow = 2, ncol = 2, byrow = TRUE)
LowAdj <- 0.05
HighAdj <- 0.95

LineWidth <- 2
Cols <- brewer.pal(n = 3, "Dark2")
SlowCol <- Cols[1]
MedCol <- Cols[3]
FastCol <- Cols[2]

# Make the extinction graph
pdf(file = "Extinction.pdf", width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = InnerMar, oma = OuterMar)
     for(i in 1:9){
          plot(NA, NA, xlim = c(1, 200), ylim = c(0, 1), xlab = "", ylab = "",
               main = "", xaxt = "n", yaxt = "n")
          lines(x = 1:200, y = slow[SimSeq[i],], lwd = LineWidth, col = SlowCol)
          lines(x = 1:200, y = med[SimSeq[i],], lwd = LineWidth, col = MedCol)
          lines(x = 1:200, y = fast[SimSeq[i],], lwd = LineWidth, col = FastCol)
          
          # Add the axes
          axis(1, cex.axis = AxisSize)
          axis(2, cex.axis = AxisSize, las = 1)
          
          # Add the selection and environmental gradient arrows
          if(i == 1){
               arrows(x0 = TopArrow[1,1], y0 = TopArrow[1,2], x1 = TopArrow[2,1], 
                      y1 = TopArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
               arrows(x0 = SideArrow[1,1], y0 = SideArrow[1,2], x1 = SideArrow[2,1], 
                      y1 = SideArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
          }
          if(i == 3){
               legend("bottomright", bty = "n", legend = c("Slow", "Medium", "Fast"),
                      lty = 1, lwd = LineWidth, col = c(SlowCol, MedCol, FastCol))
          }
     }

     # Add the x and y axis labels
     mtext("Time", side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
     mtext("Cumulative extinction probability", side = 2, outer = TRUE, line = yLabLine, cex = TextSize)
     
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
