# This script will create the graphs from the abundance data.
SpeedIndex <- 1
SuccessIndex <- 3

SuccessWords <- c("Survived", "Extinct", "All")
source("~/Desktop/RangeShifts/ShiftingSlopesCode/ShiftingRange/ShiftParams.R")
RangeParams <- read.csv("~/Desktop/RangeShifts/ShiftingSlopesOther/RangeParameters.csv")

# Set the working directory
setwd(paste("~/Desktop/RangeShifts/ShiftingSlopesOther/ShiftingRange/", 
            SpeedWords[SpeedIndex], sep = ""))
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
LocLabels <- seq(-60, 60 + 100*SpeedNums[SpeedIndex], length.out = 6)
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
#nMessageX <- 0.8
#nMessageY <- 0.1
#MessageSize <- 1.25

# Create some useful objects for graphing the position of the range center
LengthShift <- 100
StatRangeExtent <- 121
ShiftRangeExtent <- StatRangeExtent + LengthShift*SpeedNums[SpeedIndex]
ZeroPos <- ceiling(StatRangeExtent / 2)
EndShift <- ZeroPos + LengthShift*SpeedNums[SpeedIndex]
DecimalZero <- ZeroPos / ShiftRangeExtent
DecimalEndShift <- EndShift / ShiftRangeExtent

# Create a function to convert the matrix of relative locations to absolute
#    locations
RelToAbsolute <- function(RelMat, BurnIn, LengthShift, BurnOut, v, ShiftRangeExtent,
                          BetaInit, eta){
     # Create a new matrix for the absolute locations of the appropriate size
     AbsMat <- matrix(NA, nrow = ShiftRangeExtent, ncol = ncol(RelMat))
     
     # Establish the location of beta throughout the simulation
     BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
                                eta = eta, v = v) / eta
     BetaCoord <- c(rep(BetaInit, BurnIn), BetaShift, rep(BetaShift[LengthShift], 
                                                          BurnOut))
     # Step through each generation and populate the absolute matrix
     for(g in 1:ncol(RelMat)){
          AbsXcoord <- 1:ShiftRangeExtent + BetaCoord[g]
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
load(paste(SpeedWords[SpeedIndex], "ShiftingTraitResults.rdata", sep = ""))

# Change the mean dispersal phenotype results to the log scale
SectorDisp[[SuccessIndex]][,,TimeSeq,1] <- log(SectorDisp[[SuccessIndex]][,,TimeSeq,1], base = 10)

MaxSectorFit <- rep(NA, 3)
MaxSectorDisp <- rep(NA, 3)
MinSectorFit <- rep(NA, 3)
MinSectorDisp <- rep(NA, 3)


for(i in 1:3){
     MaxSectorFit[i] <- max(SectorFit[[SuccessIndex]][,,,i], na.rm = TRUE)
     MaxSectorDisp[i] <- max(SectorDisp[[SuccessIndex]][,,,i], na.rm = TRUE)
     
     MinSectorFit[i] <- min(SectorFit[[SuccessIndex]][,,,i], na.rm = TRUE)
     MinSectorDisp[i] <- min(SectorDisp[[SuccessIndex]][,,,i], na.rm = TRUE)
}
MaxAmongFit <- max(AmongVarFit[[SuccessIndex]], na.rm = TRUE)
MaxAmongDisp <- max(AmongVarDisp[[SuccessIndex]], na.rm = TRUE)
MaxWithinFit <- max(WithinVarFit[[SuccessIndex]], na.rm = TRUE)
MaxWithinDisp <- max(WithinVarDisp[[SuccessIndex]], na.rm = TRUE)
MinAmongFit <- min(AmongVarFit[[SuccessIndex]], na.rm = TRUE)
MinAmongDisp <- min(AmongVarDisp[[SuccessIndex]], na.rm = TRUE)
MinWithinFit <- min(WithinVarFit[[SuccessIndex]], na.rm = TRUE)
MinWithinDisp <- min(WithinVarDisp[[SuccessIndex]], na.rm = TRUE)

FitSectorCols <- array(NA, dim = c(3,2,10000))
DispSectorCols <- array(NA, dim = c(3,2,10000))
FitAmongCols <- array(NA, dim = c(2,10000))
DispAmongCols <- array(NA, dim = c(2,10000))
FitWithinCols <- array(NA, dim = c(2,10000))
DispWithinCols <- array(NA, dim = c(2,10000))

for(i in 1:3){
     FitSectorCols[i,1,] <- seq(MinSectorFit[i], MaxSectorFit[i], length.out = 10000)
     DispSectorCols[i,1,] <- seq(MinSectorDisp[i], MaxSectorDisp[i], length.out = 10000)
     
     FitSectorCols[i,2,] <- jet.col(10000)
     DispSectorCols[i,2,] <- jet.col(10000)
}
FitAmongCols[1,] <- seq(MinAmongFit, MaxAmongFit, length.out = 10000)
DispAmongCols[1,] <- seq(MinAmongDisp, MaxAmongDisp, length.out = 10000)
FitWithinCols[1,] <- seq(MinWithinFit, MaxWithinFit, length.out = 10000)
DispWithinCols[1,] <- seq(MinWithinDisp, MaxWithinDisp, length.out = 10000)
FitAmongCols[2,] <- jet.col(10000)
DispAmongCols[2,] <- jet.col(10000)
FitWithinCols[2,] <- jet.col(10000)
DispWithinCols[2,] <- jet.col(10000)

#if(SuccessIndex == 1){
#     NumSims <- rowSums(Success)
#} else if(SuccessIndex == 2){
#     NumSims <- ncol(Sucess) - rowSums(Success)
#} else if(SuccessIndex == 3){
#     NumSims <- ncol(Success)
#}

# --------------------- Now use a loop to make a mean, among replicate variance,
#    and within replicate variance graph for both fitness and dispersal for each
#    summary value (expected trait value, genetic variance, phenotypic variance)
PlotNames <- c("Mu", "GenVar", "PhenVar")
for(v in 1:3){
     ############# First the fitness graphs
     Fit <- c(paste("Fitness/", SuccessWords[SuccessIndex], "FitSector", PlotNames[v], ".pdf", sep = ""),
              paste("Fitness/", SuccessWords[SuccessIndex], "FitAmong", PlotNames[v], ".pdf", sep = ""),
              paste("Fitness/", SuccessWords[SuccessIndex], "FitWithin", PlotNames[v], ".pdf", sep = ""))
     # Sector mean
     pdf(file = Fit[1], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               # Find the color range for the current plot and make the figure
               ColRange <- FindRange(minimum = min(SectorFit[[SuccessIndex]][i,,TimeSeq,v], na.rm = TRUE), 
                                     maximum = max(SectorFit[[SuccessIndex]][i,,TimeSeq,v], na.rm = TRUE),
                                     sequence = FitSectorCols[v,1,])
               AbsMat <- RelToAbsolute(RelMat = SectorFit[[SuccessIndex]][i,,TimeSeq,v], BurnIn = BurnIn,
                                       LengthShift = LengthShift, BurnOut = BurnOut, v = SpeedNums[SpeedIndex],
                                       ShiftRangeExtent = ShiftRangeExtent, BetaInit = BetaInit, 
                                       eta = RangeParams$eta[1])
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
               
               # Add in the number of simulations
               #SimMessage <- paste("n = ", NumSims[i], sep = "")
               #text(x = nMessageX, y = nMessageY, labels = SimMessage, cex = MessageSize)
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
                    ColRange <- FindRange(minimum = min(AmongVarFit[[SuccessIndex]][i,,TimeSeq], na.rm = TRUE), 
                                          maximum = max(AmongVarFit[[SuccessIndex]][i,,TimeSeq], na.rm = TRUE),
                                          sequence = FitAmongCols[1,])
                    AbsMat <- RelToAbsolute(RelMat = AmongVarFit[[SuccessIndex]][i,,TimeSeq], BurnIn = BurnIn,
                                            LengthShift = LengthShift, BurnOut = BurnOut, v = SpeedNums[SpeedIndex],
                                            ShiftRangeExtent = ShiftRangeExtent, BetaInit = BetaInit, 
                                            eta = RangeParams$eta[1])
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
                    
                    # Add in the number of simulations
                    #SimMessage <- paste("n = ", NumSims[i], sep = "")
                    #text(x = nMessageX, y = nMessageY, labels = SimMessage, cex = MessageSize)
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
                    ColRange <- FindRange(minimum = min(WithinVarFit[[SuccessIndex]][i,,TimeSeq], na.rm = TRUE), 
                                          maximum = max(WithinVarFit[[SuccessIndex]][i,,TimeSeq], na.rm = TRUE),
                                          sequence = FitWithinCols[1,])
                    AbsMat <- RelToAbsolute(RelMat = WithinVarFit[[SuccessIndex]][i,,TimeSeq], BurnIn = BurnIn,
                                            LengthShift = LengthShift, BurnOut = BurnOut, v = SpeedNums[SpeedIndex],
                                            ShiftRangeExtent = ShiftRangeExtent, BetaInit = BetaInit, 
                                            eta = RangeParams$eta[1])
                    image2D(z = AbsMat, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                            main = "", col = FitWithinCols[2,ColRange], colkey = FALSE)
          
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
                    
                    # Add in the number of simulations
                    #SimMessage <- paste("n = ", NumSims[i], sep = "")
                    #text(x = nMessageX, y = nMessageY, labels = SimMessage, cex = MessageSize)
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
     Disp <- c(paste("Dispersal/", SuccessWords[SuccessIndex], "DispSector", PlotNames[v], ".pdf", sep = ""), 
              paste("Dispersal/", SuccessWords[SuccessIndex], "DispAmong", PlotNames[v], ".pdf", sep = ""),   
              paste("Dispersal/", SuccessWords[SuccessIndex], "DispWithin", PlotNames[v], ".pdf", sep = ""))  
     # Sector mean
     pdf(file = Disp[1], width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               # Find the color range for the current plot and make the figure
               ColRange <- FindRange(minimum = min(SectorDisp[[SuccessIndex]][i,,TimeSeq,v], na.rm = TRUE), 
                                     maximum = max(SectorDisp[[SuccessIndex]][i,,TimeSeq,v], na.rm = TRUE),
                                     sequence = DispSectorCols[v,1,])
               AbsMat <- RelToAbsolute(RelMat = SectorDisp[[SuccessIndex]][i,,TimeSeq,v], BurnIn = BurnIn,
                                       LengthShift = LengthShift, BurnOut = BurnOut, v = SpeedNums[SpeedIndex],
                                       ShiftRangeExtent = ShiftRangeExtent, BetaInit = BetaInit, 
                                       eta = RangeParams$eta[1])
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
               
               # Add in the number of simulations
               #SimMessage <- paste("n = ", NumSims[i], sep = "")
               #text(x = nMessageX, y = nMessageY, labels = SimMessage, cex = MessageSize)
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
                    ColRange <- FindRange(minimum = min(AmongVarDisp[[SuccessIndex]][i,,TimeSeq], na.rm = TRUE), 
                                          maximum = max(AmongVarDisp[[SuccessIndex]][i,,TimeSeq], na.rm = TRUE),
                                          sequence = DispAmongCols[1,])
                    AbsMat <- RelToAbsolute(RelMat = AmongVarDisp[[SuccessIndex]][i,,TimeSeq], BurnIn = BurnIn,
                                            LengthShift = LengthShift, BurnOut = BurnOut, v = SpeedNums[SpeedIndex],
                                            ShiftRangeExtent = ShiftRangeExtent, BetaInit = BetaInit, 
                                            eta = RangeParams$eta[1])
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
                    
                    # Add in the number of simulations
                    #SimMessage <- paste("n = ", NumSims[i], sep = "")
                    #text(x = nMessageX, y = nMessageY, labels = SimMessage, cex = MessageSize)
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
                    ColRange <- FindRange(minimum = min(WithinVarDisp[[SuccessIndex]][i,,TimeSeq], na.rm = TRUE), 
                                          maximum = max(WithinVarDisp[[SuccessIndex]][i,,TimeSeq], na.rm = TRUE),
                                          sequence = DispWithinCols[1,])
                    AbsMat <- RelToAbsolute(RelMat = WithinVarDisp[[SuccessIndex]][i,,TimeSeq], BurnIn = BurnIn,
                                            LengthShift = LengthShift, BurnOut = BurnOut, v = SpeedNums[SpeedIndex],
                                            ShiftRangeExtent = ShiftRangeExtent, BetaInit = BetaInit, 
                                            eta = RangeParams$eta[1])
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
                    
                    # Add in the number of simulations
                    #SimMessage <- paste("n = ", NumSims[i], sep = "")
                    #text(x = nMessageX, y = nMessageY, labels = SimMessage, cex = MessageSize)
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
