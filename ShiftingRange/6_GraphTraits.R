# This script will create the graphs from the abundance data.

# Set the working directory
setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/StationaryRange/")
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
LocLabels <- seq(-75, 75, by = 30)
TimeSeq <- 1:2000
TimeLabels <- seq(0, 2000, by = 400)
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
load("StationaryTraitResults.rdata")
MaxSectorFit <- rep(NA, 3)
MaxSectorDisp <- rep(NA, 3)
MaxAmongFit <- rep(NA, 3)
MaxAmongDisp <- rep(NA, 3)
MaxWithinFit <- rep(NA, 3)
MaxAmongFit <- rep(NA, 3)
for(i in 1:3){
     MaxSectorFit[i] <- max(SectorFit[,,,i])
     MaxSectorDisp[i] <- max(SectorDisp[,,,i])
     MaxAmongFit[i] <- max(AmongVarFit[,,,i])
     MaxAmongDisp[i] <- max(AmongVarDisp[,,,i])
     MaxWithinFit[i] <- max(WithinVarFit[,,,i])
     MaxWithinDisp[i] <- max(WithinVarDisp[,,,i])
}

FitSectorCols <- array(NA, dim = c(3,2,10000))
DispSectorCols <- array(NA, dim = c(3,2,10000))
FitAmongCols <- array(NA, dim = c(3,2,10000))
DispAmongCols <- array(NA, dim = c(3,2,10000))
FitWithinCols <- array(NA, dim = c(3,2,10000))
DispWithinCols <- array(NA, dim = c(3,2,10000))

for(i in 1:3){
     FitSectorCols[i,1,] <- seq(0, MaxSectorFit[i], length.out = 10000)
     DispSectorCols[i,1,] <- seq(0, MaxSectorDisp[i], length.out = 10000)
     FitAmongCols[i,1,] <- seq(0, MaxAmongFit[i], length.out = 10000)
     DispAmongCols[i,1,] <- seq(0, MaxAmongDisp[i], length.out = 10000)
     FitWithinCols[i,1,] <- seq(0, MaxWithinFit[i], length.out = 10000)
     DispWithinCols[i,1,] <- seq(0, MaxWithinDisp[i], length.out = 10000)
     
     FitSectorCols[i,2,] <- jet.col(10000)
     DispSectorCols[i,2,] <- jet.col(10000)
     FitAmongCols[i,2,] <- jet.col(10000)
     DispAmongCols[i,2,] <- jet.col(10000)
     FitWithinCols[i,2,] <- jet.col(10000)
     DispWithinCols[i,2,] <- jet.col(10000)
}

# --------------------- Now use a loop to make a mean, among replicate variance,
#    and within replicate variance graph for both fitness and dispersal for each
#    summary value (expected trait value, genetic variance, phenotypic variance)
PlotNames <- c("Mu", "GenVar", "PhenVar")
for(v in 1:3){
     ############# First the fitness graphs
     Fit <- c(paste("FitSector", PlotNames[v], ".pdf", sep = ""),
              paste("FitAmong", PlotNames[v], ".pdf", sep = ""),
              paste("FitWithin", PlotNames[v], ".pdf", sep = ""))
     # Sector mean
     pdf(file = Fit[1], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               # Find the color range for the current plot and make the figure
               ColRange <- FindRange(minimum = min(SectorFit[i,,TimeSeq,v]), 
                                     maximum = max(SectorFit[i,,TimeSeq,v]),
                                     sequence = FitSectorCols[v,1,])
               image2D(z = SectorFit[i,,TimeSeq,v], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                       main = "", col = FitSectorCols[v,2,ColRange], colkey = FALSE)
          
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
          }
     
          # Add the color key
          colkey(col = FitSectorCols[v,2,], clim = c(0, MaxSectorFit[v]), cex.axis = AxisSize, width = ColKeyWidth)
     
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

     # Among var
     pdf(file = Fit[2], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               # Find the color range for the current plot and make the figure
               ColRange <- FindRange(minimum = min(AmongVarFit[i,,TimeSeq,v]), 
                                     maximum = max(AmongVarFit[i,,TimeSeq,v]),
                                     sequence = FitAmongCols[v,1,])
               image2D(z = AmongVarFit[i,,TimeSeq,v], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                       main = "", col = FitAmongCols[v,2,ColRange], colkey = FALSE)
          
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
          }
     
          # Add the color key
          colkey(col = FitAmongCols[v,2,], clim = c(0, MaxAmongFit[v]), cex.axis = AxisSize, width = ColKeyWidth)
     
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
     
     # Within var
     pdf(file = Fit[3], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               # Find the color range for the current plot and make the figure
               ColRange <- FindRange(minimum = min(WithinVarFit[i,,TimeSeq,v]), 
                                     maximum = max(WithinVarFit[i,,TimeSeq,v]),
                                     sequence = FitWithinCols[v,1,])
               image2D(z = WithinVarFit[i,,TimeSeq,v], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                       main = "", col = FitWithinCols[v,2,ColRange], colkey = FALSE)
          
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
          }
     
          # Add the color key
          colkey(col = FitWithinCols[v,2,], clim = c(0, MaxWithinFit[v]), cex.axis = AxisSize, width = ColKeyWidth)
     
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

     
     ############# Now the dispersal graphs
     Disp <- c(paste("DispSector", PlotNames[v], ".pdf", sep = ""),
              paste("DispAmong", PlotNames[v], ".pdf", sep = ""),
              paste("DispWithin", PlotNames[v], ".pdf", sep = ""))
     # Sector mean
     pdf(file = Disp[1], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               # Find the color range for the current plot and make the figure
               ColRange <- FindRange(minimum = min(SectorDisp[i,,TimeSeq,v]), 
                                     maximum = max(SectorDisp[i,,TimeSeq,v]),
                                     sequence = DispSectorCols[v,1,])
               image2D(z = SectorDisp[i,,TimeSeq,v], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                       main = "", col = DispSectorCols[v,2,ColRange], colkey = FALSE)
          
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
          }
     
          # Add the color key
          colkey(col = DispSectorCols[v,2,], clim = c(0, MaxSectorDisp[v]), cex.axis = AxisSize, width = ColKeyWidth)
     
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
     
     # Among var
     pdf(file = Disp[2], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               # Find the color range for the current plot and make the figure
               ColRange <- FindRange(minimum = min(AmongVarDisp[i,,TimeSeq,v]), 
                                     maximum = max(AmongVarDisp[i,,TimeSeq,v]),
                                     sequence = DispAmongCols[v,1,])
               image2D(z = AmongVarDisp[i,,TimeSeq,v], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                       main = "", col = DispAmongCols[v,2,ColRange], colkey = FALSE)
          
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
          }
     
          # Add the color key
          colkey(col = DispAmongCols[v,2,], clim = c(0, MaxAmongDisp[v]), cex.axis = AxisSize, width = ColKeyWidth)
     
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
     
     # Within var
     pdf(file = Disp[3], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               # Find the color range for the current plot and make the figure
               ColRange <- FindRange(minimum = min(WithinVarDisp[i,,TimeSeq,v]), 
                                     maximum = max(WithinVarDisp[i,,TimeSeq,v]),
                                     sequence = DispWithinCols[v,1,])
               image2D(z = WithinVarDisp[i,,TimeSeq,v], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                       main = "", col = DispWithinCols[v,2,ColRange], colkey = FALSE)
          
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
          }
     
          # Add the color key
          colkey(col = DispWithinCols[v,2,], clim = c(0, MaxWithinDisp[v]), cex.axis = AxisSize, width = ColKeyWidth)
     
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
}