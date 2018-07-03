# This script will make an appropriate graph of the initial dispersal values
#    extracted by the InitDispVals.R script.

setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/")
load("FitData.rdata")
CurSpeed <- 3
SpeedWords <- c("Slow", "Main", "Fast")

# Sort the disp values into the appropriate lists
NumSims <- 200
ExtinctSims <- vector(length = 9, mode = "list")
ExtantSims <- vector(length = 9, mode = "list")
NumExtant <- rep(0, 9)
xMin <- 0
xMax <- 0
for(p in 1:9){
     ExtinctVals <- NULL
     ExtantVals <- NULL
     for(i in 1:NumSims){
          if(FitList[[p]][[i]]$Ext[CurSpeed] == 0){
               ExtinctVals <- rbind(ExtinctVals, FitList[[p]][[i]]$w)
          } else{
               ExtantVals <- rbind(ExtantVals, FitList[[p]][[i]]$w)
          }
     }
     ExtinctResults <- data.frame(x = ExtinctVals[,1], w = ExtinctVals[,2])
     ExtantResults <- data.frame(x = ExtantVals[,1], w = ExtantVals[,2])
     ExtinctSims[[p]] <- ExtinctResults
     ExtantSims[[p]] <- ExtantResults
     xMin <- min(c(xMin, ExtinctResults$x, ExtantResults$x))
     xMax <- max(c(xMax, ExtinctResults$x, ExtantResults$x))
}
xRange <- c(xMin, xMax)

# Now make a list for both extinct and extant sims with data frames for plotting
#    NOTE: the plots will likely exclude the no local adaptation scenario, but
#    we calculate them here to confirm that
ExtantPlot <- vector(length = 9, mode = "list")
ExtinctPlot <- vector(length = 9, mode = "list")
for(p in 1:9){
     ExtantxVals <- unique(ExtantSims[[p]]$x)
     ExtinctxVals <- unique(ExtinctSims[[p]]$x)
     ExtinctPlotData <- data.frame(x = rep(NA, length(ExtinctxVals)),
                                  wBar = rep(NA, length(ExtinctxVals)),
                                  lwr = rep(NA, length(ExtinctxVals)),
                                  upr = rep(NA, length(ExtinctxVals)),
                                  n = rep(NA, length(ExtinctxVals)))
     if(!is.null(ExtantxVals)){
          ExtantPlotData <- data.frame(x = rep(NA, length(ExtantxVals)),
                                       wBar = rep(NA, length(ExtantxVals)),
                                       lwr = rep(NA, length(ExtantxVals)),
                                       upr = rep(NA, length(ExtantxVals)),
                                       n = rep(NA, length(ExtantxVals)))
          for(i in 1:length(ExtantxVals)){
               CurData <- subset(ExtantSims[[p]], x == ExtantxVals[i])
               ExtantPlotData$x[i] <- ExtantxVals[i]
               ExtantPlotData$wBar[i] <- mean(CurData$w)
               ExtantPlotData$lwr[i] <- quantile(CurData$w)[2]
               ExtantPlotData$upr[i] <- quantile(CurData$w)[4]
               ExtantPlotData$n[i] <- nrow(CurData)
          }
          ExtantPlot[[p]] <- ExtantPlotData
     }     
     for(i in 1:length(ExtinctxVals)){
          CurData <- subset(ExtinctSims[[p]], x == ExtinctxVals[i])
          ExtinctPlotData$x[i] <- ExtinctxVals[i]
          ExtinctPlotData$wBar[i] <- mean(CurData$w)
          ExtinctPlotData$lwr[i] <- quantile(CurData$w)[2]
          ExtinctPlotData$upr[i] <- quantile(CurData$w)[4]
          ExtinctPlotData$n[i] <- nrow(CurData)
     }
     ExtinctPlot[[p]] <- ExtinctPlotData
}

# Set some graphical parameters to use for the subsequent figures
yRange <- c(0, 1)
xRange <- c(-40, 40)
SmallTicks <- seq(-40, 40, by = 5)
OuterMar <- c(4, 10, 5, 2)
InnerMar <- c(1.5, 1.5, 2, 1.5)
TextSize <- 1.5
AxisSize <- 1.05
ArrowWidth <- 3
ArrowLength <- 0.25
FigWidth <- 8
FigHeight <- 4
FitAxisLabel <- expression(paste("Initial ", italic("w"), sep = ""))
SimSeq <- c(7,8,9,4,5,6) # Exclude the scenarios with no local adaptation
xLabLine <- 2
yLabLine <- 1.75
GradLabLine <- 3
GradSubLine <- 1.5
AdaptLine <- 6
TopArrow <- matrix(c(-10, 1.3, 220, 1.3), nrow = 2, ncol = 2, byrow = TRUE)
LowAdj <- 0.05
HighAdj <- 0.95
ExtinctCol <- "darkred"
ExtantCol <- "darkblue"
LineWidth <- 1

PlotName <- paste(SpeedWords[CurSpeed], "InitFitSpace.pdf", sep = "")
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     par(mfrow = c(2,3), mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          if(!is.null(ExtantPlot[[i]])){
               plot(x = ExtantPlot[[i]]$x, y = ExtantPlot[[i]]$wBar, pch = 1, col = ExtantCol,
                    lwd = LineWidth, main = "", xlab = "", ylab = "",
                    xlim = xRange, ylim = yRange, las = 1, cex.axis = AxisSize)
               points(x = ExtinctPlot[[i]]$x, y = ExtinctPlot[[i]]$wBar, pch = 1, col = ExtinctCol,
                      lwd = LineWidth)
               axis(1, at = SmallTicks, labels = FALSE, tcl = -0.25)
               # Add the error bars
               segments(x0 = ExtantPlot[[i]]$x, y0 = ExtantPlot[[i]]$lwr, x1 = ExtantPlot[[i]]$x,
                        y1 = ExtantPlot[[i]]$upr, lwd = LineWidth, col = ExtantCol)
               segments(x0 = ExtinctPlot[[i]]$x, y0 = ExtinctPlot[[i]]$lwr, x1 = ExtinctPlot[[i]]$x,
                        y1 = ExtinctPlot[[i]]$upr, lwd = LineWidth, col = ExtinctCol)
          } else{
               plot(x = ExtinctPlot[[i]]$x, y = ExtinctPlot[[i]]$wBar, pch = 1, col = ExtinctCol,
                    lwd = LineWidth, main = "", xlab = "", ylab = "",
                    xlim = xRange, ylim = yRange, las = 1, cex.axis = AxisSize)
               axis(1, at = SmallTicks, labels = FALSE, tcl = -0.25)
               segments(x0 = ExtinctPlot[[i]]$x, y0 = ExtinctPlot[[i]]$lwr, x1 = ExtinctPlot[[i]]$x,
                        y1 = ExtinctPlot[[i]]$upr, lwd = LineWidth, col = ExtinctCol)
          }
          # Add the selection and environmental gradient arrows
          if(i == 7){
               arrows(x0 = TopArrow[1,1], y0 = TopArrow[1,2], x1 = TopArrow[2,1], 
                      y1 = TopArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
               mtext("High \n adaptation", side = 2, line = AdaptLine, cex = TextSize)
          }
          if(i == 4){
               mtext("Low \n adaptation", side = 2, line = AdaptLine, cex = TextSize)
          }
          if(i == 6){
               legend("bottomright", legend = c("Extant", "Extinct"), pch = 1,
                      col = c(ExtantCol, ExtinctCol), bty = "n", cex = TextSize)
          }
     }

     # Add the x and y axis labels
     mtext("Patch position", side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
     mtext(FitAxisLabel, side = 2, outer = TRUE, line = yLabLine, cex = TextSize)

     # Add the environmental gradient text
     mtext("Gradient at range edge", side = 3, outer = TRUE, line = GradLabLine,
           cex = TextSize)
     mtext("Gradual", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
           adj = LowAdj)
     mtext("Severe", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
           adj = HighAdj)
dev.off()


