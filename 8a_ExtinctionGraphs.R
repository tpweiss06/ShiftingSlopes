# This script will create the graphs from the extinction data.
setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/")
library(RColorBrewer)
RangeParams <- read.csv("RangeParameters.csv")
NumSims <- 200

# Load the extinction data
load("SimData/Extinctions.rdata")
for(i in 1:3){
     for(p in 1:9){
          Extinctions[i,p,] <- cumsum(Extinctions[i,p,]) / NumSims
     }
}

# Set some graphical parameters to use for the subsequent figures
OuterMar <- c(1,1,0,0)
TextSize <- 1.15
LegTextSize <- 1.5
AxisSize <- 1.25
LetterPos <- c(2, 0.975)
FigWidth <- 8
FigHeight <- 3
xRange <- c(0,100)
yRange <- c(0,1)
TimeAxisSeq1 <- seq(0, 100, by = 20)
TimeAxisSea2 <- seq(0, 100, by = 5)
ExtAxisSeq1 <- seq(0, 1, by = 0.2)
ExtAxisSeq2 <- seq(0, 1, by = 0.05)
GradLineType <- c(1,2,3)
ExtAxisLabel <- "Extinction Probability"
LineWidth <- 2
Col <- "darkred"

# Make the extinction graph for all three speeds
SpeedWords <- c("Slow", "Main", "Fast")
for(i in 1:3){
     PlotName <- paste("ResultFigures/", SpeedWords[i], "Extinction.pdf", sep = "")
     pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          par(mfrow = c(1,3), oma = OuterMar)
          # No local adaptation
          plot(x = 1:100, y = Extinctions[i,1,], lty = GradLineType[1], lwd = LineWidth,
               col = Col, xaxt = "n", yaxt = "n", main = "", xlab = "", ylab = "",
               xlim = xRange, ylim = yRange, type = "l")
          lines(x = 1:100, y = Extinctions[i,2,], lty = GradLineType[2], lwd = LineWidth,
                col = Col)
          lines(x = 1:100, y = Extinctions[i,3,], lty = GradLineType[3], lwd = LineWidth,
                col = Col)
          text(x = LetterPos[1], y = LetterPos[2], labels = expression(bold("a")), cex = TextSize)
          # Add the axes
          axis(1, at = TimeAxisSeq1, cex.axis = AxisSize)
          axis(1, at = TimeAxisSea2, labels = FALSE, tcl = -0.25)
          axis(2, at = ExtAxisSeq1, las = 1, cex.axis = AxisSize)
          axis(2, at = ExtAxisSeq2, labels = FALSE, tcl = -0.25)
          # Weak local adaptation
          plot(x = 1:100, y = Extinctions[i,4,], lty = GradLineType[1], lwd = LineWidth,
               col = Col, xaxt = "n", yaxt = "n", main = "", xlab = "", ylab = "",
               xlim = xRange, ylim = yRange, type = "l")
          lines(x = 1:100, y = Extinctions[i,5,], lty = GradLineType[2], lwd = LineWidth,
                col = Col)
          lines(x = 1:100, y = Extinctions[i,6,], lty = GradLineType[3], lwd = LineWidth,
                col = Col)
          text(x = LetterPos[1], y = LetterPos[2], labels = expression(bold("b")), cex = TextSize)
          # Add the axes
          axis(1, at = TimeAxisSeq1, cex.axis = AxisSize)
          axis(1, at = TimeAxisSea2, labels = FALSE, tcl = -0.25)
          axis(2, at = ExtAxisSeq1, las = 1, cex.axis = AxisSize)
          axis(2, at = ExtAxisSeq2, labels = FALSE, tcl = -0.25)
          legend(horiz = TRUE, x = -75, y = 1.35, xpd = NA, legend = c("Shallow", "Moderate", "Stark"),
                 lty = GradLineType, col = Col, lwd = LineWidth, cex = LegTextSize, bty = "n")
          # Strong local adaptation
          plot(x = 1:100, y = Extinctions[i,7,], lty = GradLineType[1], lwd = LineWidth,
               col = Col, xaxt = "n", yaxt = "n", main = "", xlab = "", ylab = "",
               xlim = xRange, ylim = yRange, type = "l")
          lines(x = 1:100, y = Extinctions[i,8,], lty = GradLineType[2], lwd = LineWidth,
                col = Col)
          lines(x = 1:100, y = Extinctions[i,9,], lty = GradLineType[3], lwd = LineWidth,
                col = Col)
          text(x = LetterPos[1], y = LetterPos[2], labels = expression(bold("c")), cex = TextSize)
          # Add the axes
          axis(1, at = TimeAxisSeq1, cex.axis = AxisSize)
          axis(1, at = TimeAxisSea2, labels = FALSE, tcl = -0.25)
          axis(2, at = ExtAxisSeq1, las = 1, cex.axis = AxisSize)
          axis(2, at = ExtAxisSeq2, labels = FALSE, tcl = -0.25)

          # Add the x and y axis labels
          mtext("Generation", side = 1, outer = TRUE, line = -1, cex = TextSize)
          mtext(ExtAxisLabel, side = 2, outer = TRUE, line = -1, cex = TextSize)
     dev.off()
}
