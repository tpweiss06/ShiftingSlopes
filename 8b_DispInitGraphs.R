# This script will make an appropriate graph of the initial dispersal values
setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/")
load("SimData/InitDispData.rdata")
CurSpeed <- 3
SpeedWords <- c("Slow", "Main", "Fast")

# Calculate the dispersal phenotype corresponding to the current speed of climate
#    change
r <- log(2)
patch_speeds <- c(0.5, 1, 2)
eta <- 50
v <- patch_speeds[CurSpeed]*eta
DispPhen <- sqrt(v^2 / r)

# Sort the disp values into the appropriate lists
AllSims <- vector(length = 9, mode = "list")
ExtantSims <- vector(length = 9, mode = "list")
NumExtant <- rep(0, 9)
NumSims <- 200
xMin <- 0
xMax <- 0
for(p in 1:9){
     AllVals <- NULL
     ExtantVals <- NULL
     for(i in 1:NumSims){
          AllVals <- c(AllVals, DispInit[[p]][[i]]$ExpDists)
          
          if(DispInit[[p]][[i]]$Ext[CurSpeed] == 1){
               NumExtant[p] <- NumExtant[p] + 1
               ExtantVals <- c(ExtantVals, DispInit[[p]][[i]]$ExpDists)
          }
     }
     AllSims[[p]] <- log(AllVals, base = 10)
     if(!is.null(ExtantVals)){
          ExtantSims[[p]] <- log(ExtantVals, base = 10)
     }
     xMin <- min(c(xMin, AllSims[[p]]))
     xMax <- max(c(xMax, AllSims[[p]]))
}
xRange <- c(xMin, xMax)
ExtantFile <- paste("SimData/", SpeedWords[CurSpeed], "NumExtant.rdata", sep = "")
save(NumExtant, file = ExtantFile)

# Set some graphical parameters to use for the subsequent figures
options(scipen = -10)
FigMat <- matrix(NA, nrow = 3, ncol = 12)
FigMat[1,] <- c(rep(1,4), rep(2,4), rep(3,4))
FigMat[2,] <- c(rep(4,4), rep(5,4), rep(6,4))
FigMat[3,] <- c(rep(7,4), rep(8,4), rep(9,4))
OuterMar <- c(4, 12, 6, 2)
InnerMar <- c(1.5, 2.25, 1.5, 2.25)
TextSize <- 1.5
AxisSize <- 1.15
ArrowWidth <- 3
ArrowLength <- 0.25
FigWidth <- 8
FigHeight <- 6
DispAxisLabel <- expression(paste("Equilibrium ", italic("log"), "(", italic("d"), ")", sep = ""))
SimSeq <- c(7,8,9,4,5,6,1,2,3)
xLabLine <- 2
yLabLine <- 3
GradLabLine <- 3.5
GradSubLine <- 1.5
AdaptLabLine <- 9
AdaptSubLine <- 7
TopArrow <- matrix(c(-0.75, 240000, 17, 240000), nrow = 2, ncol = 2, byrow = TRUE)
SideArrow <- matrix(c(-5, -550000, -5, 175000), nrow = 2, ncol = 2, byrow = TRUE)
LowAdj <- 0.05
HighAdj <- 0.95
AllCol <- "lightskyblue"
ExtantCol <- "navyblue"
HistBreaks <- seq(-9.4, 3.4, by = 0.2)
xticks <- seq(-1, 3, by = 0.2)
PlotName <- paste("ResultFigures/", SpeedWords[CurSpeed], "InitDispVals.pdf", sep = "")
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          AllHist <- hist(AllSims[[i]], plot = FALSE, breaks = HistBreaks)
          plot(AllHist, col = AllCol, main = "", xlab = "", ylab = "", cex.axis = AxisSize,
               las = 1, xlim = xRange, xaxt = "n")
          if(NumExtant[i] > 0){
               hist(ExtantSims[[i]], col = ExtantCol, breaks = HistBreaks, add = TRUE)
          }
          # Add a vertical line at the dispersal phenotype tracking climate change
          abline(v = log(DispPhen, base = 10), lty = 2, col = "red", lwd = 2)
          options(scipen = 10)
          axis(1, cex.axis = AxisSize)
          axis(1, at = xticks, labels = FALSE, tcl = -0.25)
          # Add in the number of simulations
          SimMessage <- paste("n = ", NumExtant[i], sep = "")
          options(scipen = -10)
          if(i > 3){
               legend("topleft", legend = SimMessage, lty = 0, text.col = ExtantCol,
                      bty = "n", cex = AxisSize)
          } else{
               legend("topleft", legend = SimMessage, lty = 0, text.col = ExtantCol,
                      bty = "n", cex = AxisSize)
          }
     
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
     mtext("Gradient in niche optimum", side = 2, outer = TRUE, line = AdaptLabLine,
           cex = TextSize)
     mtext("Flat", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize,
           adj = LowAdj)
     mtext("Shallow", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize)
     mtext("Steep", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize,
           adj = HighAdj)

     mtext("Range edge", side = 3, outer = TRUE, line = GradLabLine,
           cex = TextSize)
     mtext("Gradual", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
           adj = LowAdj)
     mtext("Moderate", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize)
     mtext("Stark", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
           adj = HighAdj)
dev.off()
# Reset the option for printing integers
options(scipen = 0)
