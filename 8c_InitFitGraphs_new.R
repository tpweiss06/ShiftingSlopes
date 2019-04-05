# This script will make an appropriate graph of the initial dispersal values
setwd("~/Desktop/PostdocResearch/ShiftingSlopesOther/")
load("SimData/FitDataNew.rdata")
load("SimData/ExtSimIDs.rdata")
CurSpeed <- 3
SpeedWords <- c("Slow", "Main", "Fast")

# Use SimIDs to calculate a vector of the number of extant simulatins for
#    the current conditions
NumExtant <- rep(NA, 9)
for(i in 1:9){
     CurSims <- subset(SimIDs, Params == i)
     NumExtant[i] <- sum(CurSims[,CurSpeed+2])
}

# Set some graphical parameters to use for the subsequent figures
yRange <- c(0, 1)
xRange <- c(-50, 50)
SmallTicks <- seq(-35, 35, by = 5)
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
AdaptLabLine <- 8
AdaptSubLine <- 5
TopArrow <- matrix(c(-12, 1.3, 250, 1.3), nrow = 2, ncol = 2, byrow = TRUE)
LowAdj <- 0.15
HighAdj <- 0.85
ExtinctCol <- "red"
ExtantCol <- "navyblue"
LineWidth <- 1

PlotName <- paste("ResultFigures/", SpeedWords[CurSpeed], "InitFitSpace.pdf", sep = "")
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     par(mfrow = c(2,3), mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          if(NumExtant[i] > 0){
               plot(x = FitData[[i]][[CurSpeed]]$Extant$x, y = FitData[[i]][[CurSpeed]]$Extant$wBar,
                    pch = 1, col = ExtantCol, lwd = LineWidth, main = "", xlab = "", 
                    ylab = "", xlim = xRange, ylim = yRange, las = 1, cex.axis = AxisSize)
               points(x = FitData[[i]][[CurSpeed]]$Extinct$x, y = FitData[[i]][[CurSpeed]]$Extinct$wBar, 
                      pch = 6, col = ExtinctCol, lwd = LineWidth)
               axis(1, at = SmallTicks, labels = FALSE, tcl = -0.25)
               # Add the error bars
               segments(x0 = FitData[[i]][[CurSpeed]]$Extant$x, y0 = FitData[[i]][[CurSpeed]]$Extant$lwr, 
                        x1 = FitData[[i]][[CurSpeed]]$Extant$x, y1 = FitData[[i]][[CurSpeed]]$Extant$upr, 
                        lwd = LineWidth, col = ExtantCol)
               segments(x0 = FitData[[i]][[CurSpeed]]$Extinct$x, y0 = FitData[[i]][[CurSpeed]]$Extinct$lwr, 
                        x1 = FitData[[i]][[CurSpeed]]$Extinct$x, y1 = FitData[[i]][[CurSpeed]]$Extinct$upr, 
                        lwd = LineWidth, col = ExtinctCol)
          } else{
               plot(x = FitData[[i]][[CurSpeed]]$Extinct$x, y = FitData[[i]][[CurSpeed]]$Extinct$wBar,
                    pch = 1, col = ExtinctCol, lwd = LineWidth, main = "", xlab = "", 
                    ylab = "", xlim = xRange, ylim = yRange, las = 1, cex.axis = AxisSize)
               axis(1, at = SmallTicks, labels = FALSE, tcl = -0.25)
               segments(x0 = FitData[[i]][[CurSpeed]]$Extinct$x, y0 = FitData[[i]][[CurSpeed]]$Extinct$lwr, 
                        x1 = FitData[[i]][[CurSpeed]]$Extinct$x, y1 = FitData[[i]][[CurSpeed]]$Extinct$upr, 
                        lwd = LineWidth, col = ExtinctCol)
          }
          # Add the selection and environmental gradient arrows
          mtext("Gradient in niche optimum", side = 2, outer = TRUE, line = AdaptLabLine,
                cex = TextSize)
          mtext("Shallow", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize,
                adj = LowAdj)
          mtext("Steep", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize,
                adj = HighAdj)
          if(i == 7){
               arrows(x0 = TopArrow[1,1], y0 = TopArrow[1,2], x1 = TopArrow[2,1], 
                      y1 = TopArrow[2,2], length = ArrowLength, lwd = ArrowWidth, 
                      xpd = NA)
          }
          if(i == 6){
               legend("bottomright", legend = c("Extant", "Extinct"), pch = c(1,6),
                      col = c(ExtantCol, ExtinctCol), bty = "n", cex = TextSize)
          }
     }

     # Add the x and y axis labels
     mtext("Patch position", side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
     mtext(FitAxisLabel, side = 2, outer = TRUE, line = yLabLine, cex = TextSize)

     # Add the environmental gradient text
     mtext("Range edge", side = 3, outer = TRUE, line = GradLabLine,
           cex = TextSize)
     mtext("Gradual", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
           adj = LowAdj)
     mtext("Stark", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
           adj = HighAdj)
dev.off()


