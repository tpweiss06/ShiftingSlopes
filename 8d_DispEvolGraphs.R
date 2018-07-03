# This script will create the dispersal evolution graphs

setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/")
load("DispEvolData.rdata")

# Add a new column to the DispEvol data frame for a color variable for distance. 
#    This variable will scale from black (Dist = 100) to white (Dist = 0). Also,
#    make a variable for an extant variable (to see if it makes sense)
DispEvol$ExtantCol <- ifelse(DispEvol$Dist == max(DispEvol$Dist), "darkblue", "darkred")

# Calculate the Dispersal range for the graph axes
DispRange <- range(c(DispEvol$InitDisp, DispEvol$FinalDisp))
LogRange <- c(-1.5, 3.5)

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
DispLabels <- seq(-1.5, 3.5, length.out = 6)
DispMinorTicks <- seq(-1.5, 3.5, length.out = 21)
SimSeq <- c(7,8,9,4,5,6,1,2,3)
xLabLine <- 1.75
yLabLine <- 1.5
GradLabLine <- 3
GradSubLine <- 1.5
AdaptLabLine <- 7
AdaptSubLine <- 5
TopArrow <- matrix(c(-0.5, 4.75, 15, 4.75), nrow = 2, ncol = 2, byrow = TRUE)
SideArrow <- matrix(c(-3.75, -15, -3.75, 2), nrow = 2, ncol = 2, byrow = TRUE)
LowAdj <- 0.05
HighAdj <- 0.95

# Make the graph
PlotName <- "DispEvol.pdf"
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          ParamData <- subset(DispEvol, Param == i)
          plot(NA, NA, xlim = LogRange, ylim = LogRange, main = "", xlab = "",
               ylab = "", xaxt = "n", yaxt = "n")
          points(x = log(ParamData$InitDisp, base = 10), y = log(ParamData$FinalDisp, base = 10), 
                 pch = 21, col = ParamData$ExtantCol, bg = ParamData$Distcol)
          abline(a = 0, b = 1, lty = 2, lwd = 2, col = "grey")
          
          # Add the axes
          axis(1, at = DispLabels, cex.axis = AxisSize)
          axis(1, at = DispMinorTicks, labels = FALSE, tcl = -0.25)
          axis(2, at = DispLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = DispMinorTicks, labels = FALSE, tcl = -0.25)
          
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
     mtext("Initial dispersal", side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
     mtext("Final dispersal", side = 2, outer = TRUE, line = yLabLine, cex = TextSize)
     
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


DispEvol$Delta <- log(DispEvol$FinalDisp, base = 10) - log(DispEvol$InitDisp, base = 10)
xRange <- c(-0.3, 0.3)
AllCol <- "grey50"
ExtantCol <- "springgreen4"
DispAxisLabel <- expression(paste(Delta, "log(", bar(d), ")"), sep = "")
TopArrow <- matrix(c(-0.2, 5000, 1.75, 5000), nrow = 2, ncol = 2, byrow = TRUE)
SideArrow <- matrix(c(-0.585, -11500, -0.585, 3000), nrow = 2, ncol = 2, byrow = TRUE)
yLabLine <- 1.75
xLabLine <- 2
xTicks <- seq(-0.3, 0.3, by = 0.1)
xticks <- seq(-0.3, 0.3, by = 0.02)
xLabels <- c("-0.3", "", "-0.1", "", "0.1", "", "0.3")
HistBreaks <- seq(-0.3, 0.3, by = 0.02)
PlotName <- "DispEvolHists.pdf"
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          ParamData <- subset(DispEvol, Param == i)
          SuccessData <- subset(ParamData, Dist == 100)
          AllHist <- hist(ParamData$Delta, plot = FALSE, breaks = HistBreaks)
          plot(AllHist, col = AllCol, main = "", xlab = "", ylab = "", cex.axis = AxisSize,
               las = 1, xlim = xRange, xaxt = "n")
          axis(1, at = xTicks, labels = xLabels, cex = AxisSize)
          axis(1, at = xticks, labels = FALSE, tcl = -0.25)
          hist(SuccessData$Delta, col = ExtantCol, breaks = HistBreaks, add = TRUE)
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
     mtext(DispAxisLabel, side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
     mtext("Frequency", side = 2, outer = TRUE, line = yLabLine, cex = TextSize)

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


# Calculate the observed differences between extant and extinct simulations
#    in terms of mean, standard deviation, and skew
library(moments)
ObsMeanDiff <- rep(NA, 9)
ObsSigmaDiff <- rep(NA, 9)
ObsSkewDiff <- rep(NA, 9)
Mean_p <- rep(NA, 9)
Sigma_p <- rep(NA, 9)
Skew_p <- rep(NA, 9)
for(p in 1:9){
     ExtantData <- subset(DispEvol, (Dist == 100) & (Param == p))
     ExtinctData <- subset(DispEvol, (Dist != 100) & (Param == p))
     ObsMeanDiff[p] <- mean(ExtantData$Delta) - mean(ExtinctData$Delta)
     ObsSigmaDiff[p] <- sd(ExtantData$Delta) - sd(ExtinctData$Delta)
     ObsSkewDiff[p] <- skewness(ExtantData$Delta) - skewness(ExtinctData$Delta)
     
     N <- nrow(ExtantData) + nrow(ExtinctData)
     ParamData <- subset(DispEvol, Param == p)
     MeanDiffs <- rep(NA, 1000)
     SigmaDiffs <- rep(NA, 1000)
     SkewDiffs <- rep(NA, 1000)
     for(i in 1:1000){
          ExtantSeq <- sample(1:N, size = nrow(ExtantData), replace = FALSE)
          ExtinctSeq <- setdiff(1:N, ExtantSeq)
          ExtantPerm <- ParamData$Delta[ExtantSeq]
          ExtinctPerm <- ParamData$Delta[ExtinctSeq]
          MeanDiffs[i] <- abs(mean(ExtantPerm) - mean(ExtinctPerm))
          SigmaDiffs[i] <- abs(sd(ExtantPerm) - sd(ExtinctPerm))
          SkewDiffs[i] <- abs(skewness(ExtantPerm) - skewness(ExtinctPerm))
     }
     Mean_p[p] <- sum(MeanDiffs >= abs(ObsMeanDiff[p])) / 1000
     Sigma_p[p] <- sum(SigmaDiffs >= abs(ObsSigmaDiff[p])) / 1000
     Skew_p[p] <- sum(SkewDiffs >= abs(ObsSkewDiff[p])) / 1000
     print(p)
}
ObsMeanDiff
Mean_p
ObsSigmaDiff
Sigma_p
ObsSkewDiff
Skew_p

# Now make the spatial graph
DispEvol$Delta <- log(DispEvol$FinalDisp, base = 10) - log(DispEvol$InitDisp, base = 10)
DispAxisLabel <- expression(paste(Delta, "log(", bar(d), ")", sep = ""))
#DispEvol$Delta <- DispEvol$FinalDisp - DispEvol$InitDisp
#DispAxisLabel <- expression(paste(Delta, bar(d), sep = ""))
# Make a master list with the plot data
PlotData <- vector(mode = "list", length = 9)
yMin <- 0
yMax <- 0
for(p in 1:9){
     PlotData[[p]] <- vector(mode = "list", length = 2)
     ExtantData <- subset(DispEvol, (Param == p) & (Dist == 100))
     ExtinctData <- subset(DispEvol, (Param == p) & (Dist != 100))
     Extant_xVals <- unique(ExtantData$Relx)
     Extinct_xVals <- unique(ExtinctData$Relx)
     ExtantPlot <- data.frame(x = Extant_xVals, DeltaBar = rep(NA, length(Extant_xVals)),
                         lwr = rep(NA, length(Extant_xVals)), upr = rep(NA, length(Extant_xVals)))
     ExtinctPlot <- data.frame(x = Extinct_xVals, DeltaBar = rep(NA, length(Extinct_xVals)),
                               lwr = rep(NA, length(Extinct_xVals)), upr = rep(NA, length(Extinct_xVals)))
     for(i in 1:length(Extant_xVals)){
          xData <- subset(ExtantData, Relx == Extant_xVals[i])
          ExtantPlot$DeltaBar[i] <- mean(xData$Delta)
          ExtantPlot$lwr[i] <- quantile(xData$Delta)[2]
          ExtantPlot$upr[i] <- quantile(xData$Delta)[4]
     }
     for(i in 1:length(Extinct_xVals)){
          xData <- subset(ExtinctData, Relx == Extinct_xVals[i])
          ExtinctPlot$DeltaBar[i] <- mean(xData$Delta)
          ExtinctPlot$lwr[i] <- quantile(xData$Delta)[2]
          ExtinctPlot$upr[i] <- quantile(xData$Delta)[4]
     }
     PlotData[[p]][[1]] <- ExtantPlot
     PlotData[[p]][[2]] <- ExtinctPlot
     yMin <- min(c(yMin, min(ExtantPlot$lwr), min(ExtinctPlot$lwr)))
     yMax <- max(c(yMax, max(ExtantPlot$upr), max(ExtinctPlot$upr)))
}

xRange <- range(DispEvol$Relx)
yRange <- c(yMin, yMax)
yRangeNice <- c(-0.1, 0.06)
yAxisLabels <- seq(yRangeNice[1], yRangeNice[2], length.out = 6)
yAxisMinorTicks <- seq(yRangeNice[1], yRangeNice[2], length.out = 21)
ExtinctCol <- "darkred"
ExtantCol <- "darkblue"
TopArrow <- matrix(c(-20, 0.1, 200, 0.1), nrow = 2, ncol = 2, byrow = TRUE)
SideArrow <- matrix(c(-73, -0.6, -73, 0.02), nrow = 2, ncol = 2, byrow = TRUE)
InnerMar <- c(2, 2, 2, 2)
PlotName <- "DispEvolSpatial.pdf"
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          
          plot(PlotData[[i]][[1]]$x, PlotData[[i]][[1]]$DeltaBar, xlim = xRange, 
               ylim = yRangeNice, main = "", xlab = "", ylab = "", yaxt = "n",
               pch = 1, col = ExtantCol)
          points(PlotData[[i]][[2]]$x, PlotData[[i]][[2]]$DeltaBar, pch = 1, 
                 col = ExtinctCol)
          abline(h = 0, lty = 2, col = "grey")
          # Add the error bars
          segments(x0 = PlotData[[i]][[1]]$x, y0 = PlotData[[i]][[1]]$lwr, 
                   x1 = PlotData[[i]][[1]]$x, y1 = PlotData[[i]][[1]]$upr, 
                   col = ExtantCol)
          segments(x0 = PlotData[[i]][[2]]$x, y0 = PlotData[[p]][[2]]$lwr, 
                   x1 = PlotData[[i]][[2]]$x, y1 = PlotData[[p]][[2]]$upr, 
                   col = ExtinctCol)
          # Add the y axis
          axis(2, at = yAxisLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = yAxisMinorTicks, labels = FALSE, tcl = -0.25)
          
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
     mtext("Relative patch position", side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
     mtext(DispAxisLabel, side = 2, outer = TRUE, line = yLabLine, cex = TextSize)

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


