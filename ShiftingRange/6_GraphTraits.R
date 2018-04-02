# This script will create the graphs from the abundance data.

SpeedWord <- "Slow"
SpeedNum <- 1

# Set the working directory
setwd(paste("~/Desktop/RangeShifts/ShiftingSlopesOther/ShiftingRange/", 
            SpeedWord, sep = ""))
library(plot3D)
source("~/Desktop/RangeShifts/ShiftingSlopesCode/SimFunctions.R")

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

# Create a function to convert the matrix of relative locations to absolute
#    locations
RelToAbsolute <- function(RelMat, BurnIn, LengthShift, BurnOut, v, RangeExtent,
                          BetaInit, eta){
     # Create a new matrix for the absolute locations of the appropriate size
     AbsMat <- matrix(NA, nrow = RangeExtent + v*LengthShift, ncol = ncol(RelMat))
     
     # Establish the location of beta throughout the simulation
     BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
                                eta = eta, v = v) / 50
     BetaCoord <- c(rep(BetaInit, BurnIn), BetaShift, rep(BetaShift[LengthShift], 
                                                          BurnOut))
     # Step through each generation and populate the absolute matrix
     for(g in 1:ncol(RelMat)){
          AbsXcoord <- 1:RangeExtent + BetaCoord[g]
          vals <- which(!is.na(RelMat[,g]))
          for(i in vals){
               AbsMat[AbsXcoord[i],g] <- RelMat[i,g]
          }
     }
     return(AbsMat)
}

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
load(paste(SpeedWord, "ShiftingTraitResults.rdata", sep = ""))
MaxSectorFit <- rep(NA, 3)
MaxSectorDisp <- rep(NA, 3)

for(i in 1:3){
     MaxSectorFit[i] <- max(SectorFit[,,,i])
     MaxSectorDisp[i] <- max(SectorDisp[,,,i])
}
MaxAmongFit <- max(AmongVarFit)
MaxAmongDisp <- max(AmongVarDisp)
MaxWithinFit <- max(WithinVarFit)
MaxWithinDisp <- max(WithinVarDisp)

FitSectorCols <- array(NA, dim = c(3,2,10000))
DispSectorCols <- array(NA, dim = c(3,2,10000))
FitAmongCols <- array(NA, dim = c(2,10000))
DispAmongCols <- array(NA, dim = c(2,10000))
FitWithinCols <- array(NA, dim = c(2,10000))
DispWithinCols <- array(NA, dim = c(2,10000))

for(i in 1:3){
     FitSectorCols[i,1,] <- seq(0, MaxSectorFit[i], length.out = 10000)
     DispSectorCols[i,1,] <- seq(0, MaxSectorDisp[i], length.out = 10000)
     
     FitSectorCols[i,2,] <- jet.col(10000)
     DispSectorCols[i,2,] <- jet.col(10000)
}
FitAmongCols[1,] <- seq(0, MaxAmongFit, length.out = 10000)
DispAmongCols[1,] <- seq(0, MaxAmongDisp, length.out = 10000)
FitWithinCols[1,] <- seq(0, MaxWithinFit, length.out = 10000)
DispWithinCols[1,] <- seq(0, MaxWithinDisp, length.out = 10000)
FitAmongCols[2,] <- jet.col(10000)
DispAmongCols[2,] <- jet.col(10000)
FitWithinCols[2,] <- jet.col(10000)
DispWithinCols[2,] <- jet.col(10000)

# --------------------- Now use a loop to make a mean, among replicate variance,
#    and within replicate variance graph for both fitness and dispersal for each
#    summary value (expected trait value, genetic variance, phenotypic variance)
PlotNames <- c("Mu", "GenVar", "PhenVar")
for(v in 1:3){
     ############# First the fitness graphs
     Fit <- c(paste("Figures/", SpeedWord, "FitSector", PlotNames[v], ".pdf", sep = ""),
              paste("Figures/", SpeedWord, "FitAmong", PlotNames[v], ".pdf", sep = ""),
              paste("Figures/", SpeedWord, "FitWithin", PlotNames[v], ".pdf", sep = ""))
     # Sector mean
     pdf(file = Fit[1], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               # Find the color range for the current plot and make the figure
               ColRange <- FindRange(minimum = min(SectorFit[i,,TimeSeq,v], na.rm = TRUE), 
                                     maximum = max(SectorFit[i,,TimeSeq,v], na.rm = TRUE),
                                     sequence = FitSectorCols[v,1,])
               AbsMat <- RelToAbsolute(RelMat = SectorFit[i,,TimeSeq,v], BurnIn = 50,
                                       LengthShift = 100, BurnOut = 50, v = SpeedNum,
                                       RangeExtent = 121, BetaInit = 0, eta = 50)
               image2D(z = AbsMat, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
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
               
               # Add a black line for the center of the habitable range
               segments(x0 = DecimalZero, y0 = 0, x1 = DecimalZero, y1 = 0.25,
                        lwd =  EnvLineWidth)
               segments(x0 = DecimalEndShift, y0 = 0.75, x1 = DecimalEndShift, y1 = 1,
                        lwd =  EnvLineWidth)
               segments(x0 = DecimalZero, y0 = 0.25, x1 = DecimalEndShift, y1 = 0.75,
                        lwd =  EnvLineWidth)
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

     if(v == 1){
          # Among var
          pdf(file = Fit[2], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
               layout(FigMat)
               par(mar = InnerMar, oma = OuterMar)
               for(i in SimSeq){
                    # Find the color range for the current plot and make the figure
                    ColRange <- FindRange(minimum = min(AmongVarFit[i,,TimeSeq], na.rm = TRUE), 
                                          maximum = max(AmongVarFit[i,,TimeSeq], na.rm = TRUE),
                                          sequence = FitAmongCols[1,])
                    AbsMat <- RelToAbsolute(RelMat = AmongVarFit[i,,TimeSeq], BurnIn = 50,
                                            LengthShift = 100, BurnOut = 50, v = SpeedNum,
                                            RangeExtent = 121, BetaInit = 0, eta = 50)
                    image2D(z = AbsMat, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                            main = "", col = FitAmongCols[2,ColRange], colkey = FALSE)
          
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
               colkey(col = FitAmongCols[2,], clim = c(0, MaxAmongFit), cex.axis = AxisSize, width = ColKeyWidth)
     
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
                    ColRange <- FindRange(minimum = min(WithinVarFit[i,,TimeSeq], na.rm = TRUE), 
                                          maximum = max(WithinVarFit[i,,TimeSeq], na.rm = TRUE),
                                          sequence = FitWithinCols[1,])
                    AbsMat <- RelToAbsolute(RelMat = WithinVarFit[i,,TimeSeq], BurnIn = 50,
                                            LengthShift = 100, BurnOut = 50, v = SpeedNum,
                                            RangeExtent = 121, BetaInit = 0, eta = 50)
                    image2D(z = AbsMat, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
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
               
                    # Add a black line for the center of the habitable range
                    segments(x0 = DecimalZero, y0 = 0, x1 = DecimalZero, y1 = 0.25,
                             lwd =  EnvLineWidth)
                    segments(x0 = DecimalEndShift, y0 = 0.75, x1 = DecimalEndShift, y1 = 1,
                             lwd =  EnvLineWidth)
                    segments(x0 = DecimalZero, y0 = 0.25, x1 = DecimalEndShift, y1 = 0.75,
                             lwd =  EnvLineWidth)
               }
          
               # Add the color key
               colkey(col = FitWithinCols[2,], clim = c(0, MaxWithinFit), cex.axis = AxisSize, width = ColKeyWidth)
          
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
     
     ############# Now the dispersal graphs
     Disp <- c(paste("Figures/", SpeedWord, "DispSector", PlotNames[v], ".pdf", sep = ""),
              paste("Figures/", SpeedWord, "DispAmong", PlotNames[v], ".pdf", sep = ""),
              paste("Figures/", SpeedWord, "DispWithin", PlotNames[v], ".pdf", sep = ""))
     # Sector mean
     pdf(file = Disp[1], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               # Find the color range for the current plot and make the figure
               ColRange <- FindRange(minimum = min(SectorDisp[i,,TimeSeq,v], na.rm = TRUE), 
                                     maximum = max(SectorDisp[i,,TimeSeq,v], na.rm = TRUE),
                                     sequence = DispSectorCols[v,1,])
               AbsMat <- RelToAbsolute(RelMat = SectorDisp[i,,TimeSeq,v], BurnIn = 50,
                                       LengthShift = 100, BurnOut = 50, v = SpeedNum,
                                       RangeExtent = 121, BetaInit = 0, eta = 50)
               image2D(z = AbsMat, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
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
               
               # Add a black line for the center of the habitable range
               segments(x0 = DecimalZero, y0 = 0, x1 = DecimalZero, y1 = 0.25,
                        lwd =  EnvLineWidth)
               segments(x0 = DecimalEndShift, y0 = 0.75, x1 = DecimalEndShift, y1 = 1,
                        lwd =  EnvLineWidth)
               segments(x0 = DecimalZero, y0 = 0.25, x1 = DecimalEndShift, y1 = 0.75,
                        lwd =  EnvLineWidth)
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
     
     if(v == 1){
          # Among var
          pdf(file = Disp[2], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
               layout(FigMat)
               par(mar = InnerMar, oma = OuterMar)
               for(i in SimSeq){
                    # Find the color range for the current plot and make the figure
                    ColRange <- FindRange(minimum = min(AmongVarDisp[i,,TimeSeq], na.rm = TRUE), 
                                          maximum = max(AmongVarDisp[i,,TimeSeq], na.rm = TRUE),
                                          sequence = DispAmongCols[1,])
                    AbsMat <- RelToAbsolute(RelMat = AmongVarDisp[i,,TimeSeq], BurnIn = 50,
                                            LengthShift = 100, BurnOut = 50, v = SpeedNum,
                                            RangeExtent = 121, BetaInit = 0, eta = 50)
                    image2D(z = AbsMat, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                            main = "", col = DispAmongCols[2,ColRange], colkey = FALSE)
          
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
               colkey(col = DispAmongCols[2,], clim = c(0, MaxAmongDisp), cex.axis = AxisSize, width = ColKeyWidth)
     
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
                    ColRange <- FindRange(minimum = min(WithinVarDisp[i,,TimeSeq], na.rm = TRUE), 
                                          maximum = max(WithinVarDisp[i,,TimeSeq], na.rm = TRUE),
                                          sequence = DispWithinCols[1,])
                    AbsMat <- RelToAbsolute(RelMat = WithinVarDisp[i,,TimeSeq], BurnIn = 50,
                                            LengthShift = 100, BurnOut = 50, v = SpeedNum,
                                            RangeExtent = 121, BetaInit = 0, eta = 50)
                    image2D(z = AbsMat, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                            main = "", col = DispWithinCols[2,ColRange], colkey = FALSE)
          
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
               colkey(col = DispWithinCols[2,], clim = c(0, MaxWithinDisp), cex.axis = AxisSize, width = ColKeyWidth)
     
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
}