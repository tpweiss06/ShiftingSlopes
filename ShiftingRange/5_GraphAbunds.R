# This script will create the graphs from the abundance data.

SpeedWord <- "Slow"
SpeedNum <- 1

# Set the working directory
setwd(paste("~/Desktop/RangeShifts/ShiftingSlopesOther/ShiftingRange/", 
            SpeedWord, sep = ""))
library(plot3D)

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
LocLabels <- seq(-60, 60 + 100*SpeedNum, length.out = 6)
TimeSeq <- 1:200
TimeLabels <- seq(0, 200, by = 40)
SimSeq <- c(7,8,9,4,5,6,1,2,3)
ColKeyWidth <- 15
xLabLine <- 1.75
yLabLine <- 1.5
GradLabLine <- 3
GradSubLine <- 1.5
AdaptLabLine <- 7
AdaptSubLine <- 5
TopArrow <- matrix(c(0.2, 1.2, 3.5, 1.2), nrow = 2, ncol = 2, byrow = TRUE)
SideArrow <- matrix(c(-0.45, -2.5, -0.45, 0.8), nrow = 2, ncol = 2, byrow = TRUE)
LowAdj <- 0.05
HighAdj <- 0.95
EnvLineWidth <- 1

# Create some useful objects for graphing the position of the range center
RangeExtent <- 121 + 100*SpeedNum
ZeroPos <- 61
EndShift <- ZeroPos + 100*SpeedNum
DecimalZero <- ZeroPos / RangeExtent
DecimalEndShift <- EndShift / RangeExtent

# Create a function to find the indices corresponding to the subset of the a sequence
#	which is defined as the shortest possible sequence including both the given
#	minimum and maximum
FindRange <- function(minimum, maximum, sequence){
     if(class(sequence) != 'numeric'){
          sequence <- as.numeric(sequence)
     }
     LessThan <- ifelse(sequence <= minimum, 1, 0)
     MoreThan <- ifelse(sequence >= maximum, 1, 0)
     MinIndex <- sum(LessThan)
     MaxIndex <- length(sequence) - sum(MoreThan)
     return(MinIndex:MaxIndex)
}

# Get the appropriate ranges for the abundance and sigma values and set up
#    color matrices appropriately
load(paste(SpeedWord, "ShiftingAbundResults.rdata", sep = ""))
MaxAbund <- max(SectorMean, na.rm = TRUE)
MaxWithinVar <- max(WithinVar, na.rm = TRUE)
MaxAmongVar <- max(AmongVar, na.rm = TRUE)

AbundCols <- matrix(NA, nrow = 2, ncol = 10000)
WithinCols <- matrix(NA, nrow = 2, ncol = 10000)
AmongCols <- matrix(NA, nrow = 2, ncol = 10000)

AbundCols[1,] <- seq(0, MaxAbund, length.out = 10000)
WithinCols[1,] <- seq(0, MaxWithinVar, length.out = 10000)
AmongCols[1,] <- seq(0, MaxAmongVar, length.out = 10000)
AbundCols[2,] <- jet.col(10000)
WithinCols[2,] <- jet.col(10000)
AmongCols[2,] <- jet.col(10000)

# Make the mean abundance graph
PlotName <- paste("Figures/", SpeedWord, "MeanAbunds.pdf", sep = "")
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          # Find the color range for the current plot and make the figure
          ColRange <- FindRange(minimum = min(SectorMean[i,,TimeSeq]), 
                                maximum = max(SectorMean[i,,TimeSeq]),
                                sequence = AbundCols[1,])
          image2D(z = SectorMean[i,,TimeSeq], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = "", col = AbundCols[2,ColRange], colkey = FALSE)
          
          # Add the axes
          axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
          axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          
          # Add the selection and environmental gradient arrows
          if(i == 7){
               arrows(x0 = TopArrow[1,1], y0 = TopArrow[1,2], x1 = TopArrow[2,1], 
                      y1 = TopArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
               arrows(x0 = SideArrow[1,1], y0 = SideArrow[1,2], x1 = SideArrow[2,1], 
                      y1 = SideArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
          }
          
          # Add a black line for the center of the habitable range
          segments(x0 = DecimalZero, y0 = 0, x1 = DecimalZero, y1 = 0.25,
                   lwd =  EnvLineWidth)
          segments(x0 = DecimalEndShift, y0 = 0.75, x1 = DecimalEndShift, y1 = 1,
                   lwd =  EnvLineWidth)
          segments(x0 = DecimalZero, y0 = 0.25, x1 = DecimalEndShift, y1 = 0.75,
                   lwd =  EnvLineWidth)
     }

     # Add the color key
     colkey(col = AbundCols[2,], clim = c(0, MaxAbund), cex.axis = AxisSize, width = ColKeyWidth)

     # Add the x and y axis labels
     mtext("Spatial location (x)", side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
     mtext("Generation", side = 2, outer = TRUE, line = yLabLine, cex = TextSize)
     
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


# Make the within variance graph
PlotName <- paste("Figures/", SpeedWord, "AbundWithinVar.pdf", sep = "")
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          # Find the color range for the current plot and make the figure
          ColRange <- FindRange(minimum = min(WithinVar[i,,TimeSeq]), 
                                maximum = max(WithinVar[i,,TimeSeq]),
                                sequence = WithinCols[1,])
          image2D(z = WithinVar[i,,TimeSeq], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = "", col = WithinCols[2,ColRange], colkey = FALSE)
     
          # Add the axes
          axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
          axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     
          # Add the selection and environmental gradient arrows
          if(i == 7){
               arrows(x0 = TopArrow[1,1], y0 = TopArrow[1,2], x1 = TopArrow[2,1], 
                      y1 = TopArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
               arrows(x0 = SideArrow[1,1], y0 = SideArrow[1,2], x1 = SideArrow[2,1], 
                      y1 = SideArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
          }
          
          # Add a black line for the center of the habitable range
          segments(x0 = DecimalZero, y0 = 0, x1 = DecimalZero, y1 = 0.25,
                   lwd =  EnvLineWidth)
          segments(x0 = DecimalEndShift, y0 = 0.75, x1 = DecimalEndShift, y1 = 1,
                   lwd =  EnvLineWidth)
          segments(x0 = DecimalZero, y0 = 0.25, x1 = DecimalEndShift, y1 = 0.75,
                   lwd =  EnvLineWidth)
     }

     # Add the color key
     colkey(col = WithinCols[2,], clim = c(0, MaxWithinVar), cex.axis = AxisSize, width = ColKeyWidth)

     # Add the x and y axis labels
     mtext("Spatial location (x)", side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
     mtext("Generation", side = 2, outer = TRUE, line = yLabLine, cex = TextSize)

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

# Make the among variance graph
PlotName <- paste("Figures/", SpeedWord, "AbundAmongVar.pdf", sep = "")
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          # Find the color range for the current plot and make the figure
          ColRange <- FindRange(minimum = min(AmongVar[i,,TimeSeq]), 
                                maximum = max(AmongVar[i,,TimeSeq]),
                                sequence = AmongCols[1,])
          image2D(z = AmongVar[i,,TimeSeq], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = "", col = AmongCols[2,ColRange], colkey = FALSE)
     
          # Add the axes
          axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
          axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     
          # Add the selection and environmental gradient arrows
          if(i == 7){
               arrows(x0 = TopArrow[1,1], y0 = TopArrow[1,2], x1 = TopArrow[2,1], 
                      y1 = TopArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
               arrows(x0 = SideArrow[1,1], y0 = SideArrow[1,2], x1 = SideArrow[2,1], 
                      y1 = SideArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
          }
          
          # Add a black line for the center of the habitable range
          segments(x0 = DecimalZero, y0 = 0, x1 = DecimalZero, y1 = 0.25,
                   lwd =  EnvLineWidth)
          segments(x0 = DecimalEndShift, y0 = 0.75, x1 = DecimalEndShift, y1 = 1,
                   lwd =  EnvLineWidth)
          segments(x0 = DecimalZero, y0 = 0.25, x1 = DecimalEndShift, y1 = 0.75,
                   lwd =  EnvLineWidth)
     }

     # Add the color key
     colkey(col = AmongCols[2,], clim = c(0, MaxAmongVar), cex.axis = AxisSize, width = ColKeyWidth)

     # Add the x and y axis labels
     mtext("Spatial location (x)", side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
     mtext("Generation", side = 2, outer = TRUE, line = yLabLine, cex = TextSize)

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

