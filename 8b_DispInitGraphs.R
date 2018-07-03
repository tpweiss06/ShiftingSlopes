# This script will make an appropriate graph of the initial dispersal values
setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/")
load("InitDispData.rdata")
CurSpeed <- 1
SpeedWords <- c("Slow", "Main", "Fast")

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
     for(i in 1:100){
          AllVals <- c(AllVals, DispList[[p]][[i]][-1])
          if(DispList[[p]][[i]][1] == 1){
               NumExtant[p] <- NumExtant[p] + 1
               ExtantVals <- c(ExtantVals, DispList[[p]][[i]][-1])
          }
     }
     AllSims[[p]] <- log(AllVals, base = 10)
     ExtantSims[[p]] <- log(ExtantVals, base = 10)
     xMin <- min(c(xMin, AllSims[[p]]))
     xMax <- max(c(xMax, AllSims[[p]]))
}
xRange <- c(xMin, xMax)

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
DispAxisLabel <- expression(paste("Initial ", italic("log"), "(", bar(italic("d")), ")", sep = ""))
SimSeq <- c(7,8,9,4,5,6,1,2,3)
xLabLine <- 2
yLabLine <- 1.75
GradLabLine <- 3
GradSubLine <- 1.5
AdaptLabLine <- 7
AdaptSubLine <- 5
TopArrow <- matrix(c(-0.5, 2750, 13, 2750), nrow = 2, ncol = 2, byrow = TRUE)
SideArrow <- matrix(c(-3.5, -5500, -3.5, 1500), nrow = 2, ncol = 2, byrow = TRUE)
LowAdj <- 0.05
HighAdj <- 0.95
AllCol <- "grey50"
ExtantCol <- "springgreen4"
HistBreaks <- seq(-1.4, 3.2, by = 0.2)
xticks <- seq(-1, 3, by = 0.2)
PlotName <- paste(SpeedWords[CurSpeed], "InitDispVals.pdf", sep = "")
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     layout(FigMat)
     par(mar = InnerMar, oma = OuterMar)
     for(i in SimSeq){
          AllHist <- hist(AllSims[[i]], plot = FALSE, breaks = HistBreaks)
          plot(AllHist, col = AllCol, main = "", xlab = "", ylab = "", cex.axis = AxisSize,
               las = 1, xlim = xRange)
          hist(ExtantSims[[i]], col = ExtantCol, breaks = HistBreaks, add = TRUE)
          axis(1, at = xticks, labels = FALSE, tcl = -0.25)
          # Add in the number of simulations
          SimMessage <- paste("n = ", NumExtant[i], sep = "")
          if(i > 3){
               legend("topleft", legend = SimMessage, lty = 0, text.col = ExtantCol,
                      bty = "n")
          } else{
               legend("topleft", legend = SimMessage, lty = 0, text.col = ExtantCol,
                      bty = "n")
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
     mtext("Potential for local adaptation", side = 2, outer = TRUE, line = AdaptLabLine,
           cex = TextSize)
     mtext("None", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize,
           adj = LowAdj)
     mtext("High", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize,
           adj = HighAdj)

     mtext("Gradient at range edge", side = 3, outer = TRUE, line = GradLabLine,
           cex = TextSize)
     mtext("Gradual", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
           adj = LowAdj)
     mtext("Severe", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
           adj = HighAdj)
dev.off()

