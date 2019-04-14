# This script contains an example of some of the parameters and lines of code
#    to create the graphs with arrows if Allison prefers that. This isn't all
#    the necessary code, but it has everything that needs to be added back in.

yArrAdjLwr <- 0.2
yArrAdjUpr <- 3
xArrAdj <- c(1.4, 1.45, 1.5)
LowAdj <- 0.05
HighAdj <- 0.95
LowAdjFit <- 0.25
HighAdjFit <- 0.75
OptLabLine <- 10
OptWordLine <- 7.5
OptLabLineFit <- 7.5
OptWordLineFit <- 5.5
OptWords <- c("Flat", "Shallow", "Steep")
OptLabel <- "Gradient in niche optimum"

# -------------------------------------------------- 1) mean dispersal phenotype
yLabel <- expression(bar(d))
load("SimData/DispDataNew_a.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- rep(Inf, 3)
xUpr <- rep(-Inf, 3)
yLwr <- rep(Inf, 3)
yUpr <- rep(-Inf, 3)
for(p in 1:9){
     for(v in 1:3){
          CurData <- DispData[[p]][[v]]
          xLwr[v] <- min(xLwr[v], min(CurData$x))
          xUpr[v] <- max(xUpr[v], max(CurData$x))
          yLwr[v] <- min(yLwr[v], min(CurData$lwr))
          yUpr[v] <- max(yUpr[v], max(CurData$upr))
     }
}

for(v in 1:3){
     FigName <- paste("ResultFigures/", SpeedWords[v], "MeanDisp.pdf", sep = "")
     pdf(file = FigName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     par(mfrow = c(3,1), oma = OuterMar, mar = InnerMar)
     for(n in 1:3){
          if(n == 3){
               plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                    main = "", xlab = "", ylab = "", cex.axis = AxisSize, las = 1)
          }else{
               plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                    main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize, las = 1)
               axis(1, labels = FALSE)
          }
          for(e in 1:3){
               CurParam <- SimMat[n,e]
               for(i in 1:length(PlotGens)){
                    CurData <- subset(DispData[[CurParam]][[v]], g == PlotGens[i])
                    CurData <- CurData[order(CurData$x),]
                    lines(x = CurData$x, y = CurData$dBar, lty = 1,
                          col = EdgeCols[e], lwd = LineWidth)
                    lines(x = CurData$x, y = CurData$lwr, lty = 3,
                          col = EdgeCols[e])
                    lines(x = CurData$x, y = CurData$upr, lty = 3,
                          col = EdgeCols[e])
                    if(e == 1){
                         text(x = BetaPos[v,i], y = GenTextAdj*yUpr[v], 
                              labels = paste("t = ", PlotGens[i], sep = ""),
                              cex = LetterSize)
                    }
               }
          }
          # Add a letter to the current plot
          #text(x = xLetterAdj * xLwr[v], y = yLetterAdj * yUpr[v], 
          #     labels = LetterSeq[n], cex = LetterSize)
          # Add a legend on the first plot
          if(n == 1){
               LegendText <- c(EdgeWords[1], "", EdgeWords[2], "", EdgeWords[3])
               LegendCols <- c(EdgeCols[1], "white", EdgeCols[2], "white", EdgeCols[3])
               legend(x = "top", legend = LegendText, lty = c(1,0,1,0,1), col = LegendCols, 
                      bty = "n", cex = LegSize, lwd = LineWidth, horiz = TRUE,
                      inset = LegInset, xpd = NA)
          }
     }
     # Add the niche optimum words
     mtext(text = OptLabel, side = 2, outer = TRUE, line = OptLabLine,
           cex = TextSize)
     mtext(text = OptWords[1], side = 2, outer = TRUE, line = OptWordLine, 
           cex = TextSize, adj = LowAdj)
     mtext(text = OptWords[2], side = 2, outer = TRUE, line = OptWordLine, 
           cex = TextSize)
     mtext(text = OptWords[3], side = 2, outer = TRUE, line = OptWordLine, 
           cex = TextSize, adj = HighAdj)
     # Add the niche optimum arrow
     corners <- par("usr")
     arrows(x0 = xArrAdj[v]*corners[1], y0 = yArrAdjLwr*corners[4], 
            x1 = xArrAdj[v]*corners[1], y1 = yArrAdjUpr*corners[4], 
            xpd = NA, code = 2, lwd = LineWidth)
     # Finally, add the axis labels
     mtext(text = xLabel, side = 1, outer = TRUE, line = xLabLine,
           cex = TextSize)
     mtext(text = yLabel, side = 2, outer = TRUE, line = yLabLine,
           cex = TextSize, las = 1)
     dev.off()
}
rm(DispData)


