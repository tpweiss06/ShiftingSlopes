# This script will create figures showing extinction probabilities and the change 
#    in fitness, dispersal, and genetic variance in both during the range shift.

setwd("~/Desktop/PostdocResearch/ShiftingSlopesOther/")
library(RColorBrewer)

# First set some common graphical parameters to be used among all the figures
# Placement
xLetterAdj <- 1.1
yLetterAdj <- 0.95
xLabLine <- 3
GenTextAdj <- 0.9
BetaPos <- matrix(c(0,50,100,0,62.5,125,0,75,150), nrow = 3, ncol = 3, byrow = TRUE)
LegInset <- -0.2
LegInsetExt <- -0.15
# Sizing
FigWidth <- 12
FigHeight <- 9
FigPointSize <- 18
LineWidth <- 2.5
# Labels and colors
xLabel <- "Space"
EdgeWords <- c("Gradual", "Moderate", "Stark")
LetterSeq <- c("a", "b", "c")
EdgeCols <- brewer.pal(n = 3, name = "Dark2")
# Miscelaneous
SimMat <- matrix(c(7,8,9,4,5,6,1,2,3), nrow = 3, ncol = 3, byrow = TRUE)
PlotGens <- c(0, 50, 100)
InnerMar <- c(0.45,0,0.45,0)
SpeedWords <- c("Slow", "Main", "Fast")

# Now set some graphical parameters specific to the three panel dispersal graphs
yLabLine <- 5
TextSize <- 1.75
AxisSize <- 2.25
LegSize <- 2.5
LetterSize <- 2.5
OuterMar <- c(5, 9, 4, 2)

# -------------------------------------------------- 1) mean dispersal phenotype
yLabel <- "Within simulation SD in mean dispersal"
load("SimData/xVarDispData.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- rep(Inf, 3)
xUpr <- rep(-Inf, 3)
yLwr <- rep(Inf, 3)
yUpr <- rep(-Inf, 3)
for(p in 1:9){
     for(v in 1:3){
          CurData <- DispData[[p]][[v]]
          xLwr[v] <- min(xLwr[v], min(CurData$x, na.rm = TRUE))
          xUpr[v] <- max(xUpr[v], max(CurData$x, na.rm = TRUE))
          yLwr[v] <- min(yLwr[v], min(CurData$lwr, na.rm = TRUE))
          yUpr[v] <- max(yUpr[v], max(CurData$upr, na.rm = TRUE))
     }
}

for(v in 1:3){
     FigName <- paste("ResultFigures/xVar/", SpeedWords[v], "MeanDisp.pdf", sep = "")
     pdf(file = FigName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          par(mfrow = c(3,1), oma = OuterMar, mar = InnerMar)
          for(n in 1:3){
               if(n == 3){
                    plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                         main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize, las = 1)
                    axis(1, padj = 0.5, cex.axis = AxisSize)
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
               text(x = xLetterAdj * xLwr[v], y = yLetterAdj * yUpr[v], 
                    labels = LetterSeq[n], cex = LetterSize)
               # Add a legend on the first plot
               if(n == 1){
                    LegendText <- c(EdgeWords[1], "", EdgeWords[2], "", EdgeWords[3])
                    LegendCols <- c(EdgeCols[1], "white", EdgeCols[2], "white", EdgeCols[3])
                    legend(x = "top", legend = LegendText, lty = c(1,0,1,0,1), col = LegendCols, 
                           bty = "n", cex = LegSize, lwd = LineWidth, horiz = TRUE,
                           inset = LegInset, xpd = NA)
               }
          }
          # Finally, add the axis labels
          mtext(text = xLabel, side = 1, outer = TRUE, line = xLabLine,
                cex = TextSize)
          mtext(text = yLabel, side = 2, outer = TRUE, line = yLabLine,
                cex = TextSize, las = 1)
     dev.off()
}
rm(DispData)

# --------------------------------------------- 2) genetic variance in dispersal
yLabel <- "Within simulation SD in mean dispersal GenVar"
load("SimData/xVarDispGenVarData.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- rep(Inf, 3)
xUpr <- rep(-Inf, 3)
yLwr <- rep(Inf, 3)
yUpr <- rep(-Inf, 3)
for(p in 1:9){
     for(v in 1:3){
          CurData <- DispGenVarData[[p]][[v]]
          xLwr[v] <- min(xLwr[v], min(CurData$x, na.rm = TRUE))
          xUpr[v] <- max(xUpr[v], max(CurData$x, na.rm = TRUE))
          yLwr[v] <- min(yLwr[v], min(CurData$lwr, na.rm = TRUE))
          yUpr[v] <- max(yUpr[v], max(CurData$upr, na.rm = TRUE))
     }
}

for(v in 1:3){
     FigName <- paste("ResultFigures/xVar/", SpeedWords[v], "DispGenVar.pdf", sep = "")
     pdf(file = FigName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          par(mfrow = c(3,1), oma = OuterMar, mar = InnerMar)
          for(n in 1:3){
               if(n == 3){
                    plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                         main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize, las = 1)
                    axis(1, padj = 0.5, cex.axis = AxisSize)
               }else{
                    plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                         main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize, las = 1)
                    axis(1, labels = FALSE)
               }
               for(e in 1:3){
                    CurParam <- SimMat[n,e]
                    for(i in 1:length(PlotGens)){
                         CurData <- subset(DispGenVarData[[CurParam]][[v]], g == PlotGens[i])
                         CurData <- CurData[order(CurData$x),]
                         lines(x = CurData$x, y = CurData$GenVar, lty = 1,
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
               text(x = xLetterAdj * xLwr[v], y = yLetterAdj * yUpr[v], 
                    labels = LetterSeq[n], cex = LetterSize)
               # Add a legend on the first plot
               if(n == 1){
                    LegendText <- c(EdgeWords[1], "", EdgeWords[2], "", EdgeWords[3])
                    LegendCols <- c(EdgeCols[1], "white", EdgeCols[2], "white", EdgeCols[3])
                    legend(x = "top", legend = LegendText, lty = c(1,0,1,0,1), col = LegendCols, 
                           bty = "n", cex = LegSize, lwd = LineWidth, horiz = TRUE,
                           inset = LegInset, xpd = NA)
               }
          }
          # Finally, add the axis labels
          mtext(text = xLabel, side = 1, outer = TRUE, line = xLabLine,
                cex = TextSize)
          mtext(text = yLabel, side = 2, outer = TRUE, line = yLabLine,
                cex = TextSize)
     dev.off()
}
rm(DispGenVarData)

############# Now update some graphical parameters for the two panel niche plots
yLabLine <- 4
TextSize <- 2
AxisSize <- 2
LegSize <- 1.75
LetterSize <- 2
OuterMar <- c(5, 6, 4, 2)
#############

#-------------------------------------------------- 3) Graph the niche genotypes
yLabel <- "Within simulation SD in mean niche values"
load("SimData/xVarNicheData.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- rep(Inf, 3)
xUpr <- rep(-Inf, 3)
yLwr <- rep(Inf, 3)
yUpr <- rep(-Inf, 3)
for(p in 1:9){
     for(v in 1:3){
          CurData <- MismatchData[[p]][[v]]
          xLwr[v] <- min(xLwr[v], min(CurData$x, na.rm = TRUE))
          xUpr[v] <- max(xUpr[v], max(CurData$x, na.rm = TRUE))
          yLwr[v] <- min(yLwr[v], min(CurData$lwr, na.rm = TRUE))
          yUpr[v] <- max(yUpr[v], max(CurData$upr, na.rm = TRUE))
     }
}

for(v in 1:3){
     FigName <- paste("ResultFigures/xVar/", SpeedWords[v], "MeanNiche.pdf", sep = "")
     pdf(file = FigName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          par(mfrow = c(2,1), oma = OuterMar, mar = InnerMar)
          for(n in 1:2){
               if(n == 2){
                    plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                         main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize, las = 1)
                    axis(1, padj = 0.5, cex.axis = AxisSize)
               }else{
                    plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                         main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize, las = 1)
                    axis(1, labels = FALSE)
               }
               for(e in 1:3){
                    CurParam <- SimMat[n,e]
                    for(i in 1:length(PlotGens)){
                         CurData <- subset(MismatchData[[CurParam]][[v]], g == PlotGens[i])
                         CurData <- CurData[order(CurData$x),]
                         lines(x = CurData$x, y = CurData$Mismatch, lty = 1,
                               col = EdgeCols[e], lwd = LineWidth)
                         lines(x = CurData$x, y = CurData$lwr, lty = 3,
                               col = EdgeCols[e])
                         lines(x = CurData$x, y = CurData$upr, lty = 3,
                               col = EdgeCols[e])
                         # Add in the Zopt line
                         lines(x = CurData$x, y = CurData$Zopt, col = "grey")
                         if(e == 1){
                              text(x = BetaPos[v,i], y = GenTextAdj*yUpr[v], 
                                   labels = paste("t = ", PlotGens[i], sep = ""),
                                   cex = TextSize)
                         }
                    }
               }
               # Add a letter to the current plot
               text(x = xLetterAdj * xLwr[v], y = yLetterAdj * yUpr[v], 
                    labels = LetterSeq[n], cex = LetterSize)
               # Add a legend on the first plot
               if(n == 1){
                    LegendText <- c(EdgeWords[1], "", EdgeWords[2], "", EdgeWords[3])
                    LegendCols <- c(EdgeCols[1], "white", EdgeCols[2], "white", EdgeCols[3])
                    legend(x = "top", legend = LegendText, lty = c(1,0,1,0,1), col = LegendCols, 
                           bty = "n", cex = LegSize, lwd = LineWidth, horiz = TRUE,
                           inset = LegInset, xpd = NA)
               }
          }
          # Finally, add the axis labels
          mtext(text = xLabel, side = 1, outer = TRUE, line = xLabLine,
                cex = TextSize)
          mtext(text = yLabel, side = 2, outer = TRUE, line = yLabLine,
                cex = TextSize)
     dev.off()
}
rm(MismatchData)

# ----------------------------------------------- 4) genetic variance in fitness
yLabel <- "Within simulation SD in mean niche GenVar"
load("SimData/xVarFitGenVarData.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- rep(Inf, 3)
xUpr <- rep(-Inf, 3)
yLwr <- rep(Inf, 3)
yUpr <- rep(-Inf, 3)
for(p in 1:9){
     for(v in 1:3){
          CurData <- FitGenVarData[[p]][[v]]
          xLwr[v] <- min(xLwr[v], min(CurData$x, na.rm = TRUE))
          xUpr[v] <- max(xUpr[v], max(CurData$x, na.rm = TRUE))
          yLwr[v] <- min(yLwr[v], min(CurData$lwr, na.rm = TRUE))
          yUpr[v] <- max(yUpr[v], max(CurData$upr, na.rm = TRUE))
     }
}

for(v in 1:3){
     FigName <- paste("ResultFigures/xVar/", SpeedWords[v], "FitGenVar.pdf", sep = "")
     pdf(file = FigName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          par(mfrow = c(2,1), oma = OuterMar, mar = InnerMar)
          for(n in 1:2){
               if(n == 2){
                    plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                         main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize, las = 1)
                    axis(1, padj = 0.5, cex.axis = AxisSize)
               }else{
                    plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                         main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize, las = 1)
                    axis(1, labels = FALSE)
               }
               for(e in 1:3){
                    CurParam <- SimMat[n,e]
                    for(i in 1:length(PlotGens)){
                         CurData <- subset(FitGenVarData[[CurParam]][[v]], g == PlotGens[i])
                         CurData <- CurData[order(CurData$x),]
                         lines(x = CurData$x, y = CurData$GenVar, lty = 1,
                               col = EdgeCols[e], lwd = LineWidth)
                         lines(x = CurData$x, y = CurData$lwr, lty = 3,
                               col = EdgeCols[e])
                         lines(x = CurData$x, y = CurData$upr, lty = 3,
                               col = EdgeCols[e])
                         if(e == 1){
                              if((v == 3) & (i == 1) & (n == 1)){
                                   text(x = BetaPos[v,i], y = (1-GenTextAdj)*yUpr[v], 
                                        labels = paste("t = ", PlotGens[i], sep = ""),
                                        cex = TextSize)
                              }else{
                                   text(x = BetaPos[v,i], y = GenTextAdj*yUpr[v], 
                                        labels = paste("t = ", PlotGens[i], sep = ""),
                                        cex = TextSize)
                              }
                         }
                    }
               }
               # Add a letter to the current plot
               text(x = xLetterAdj * xLwr[v], y = yLetterAdj * yUpr[v], 
                    labels = LetterSeq[n], cex = LetterSize)
               # Add a legend on the first plot
               if(n == 1){
                    LegendText <- c(EdgeWords[1], "", EdgeWords[2], "", EdgeWords[3])
                    LegendCols <- c(EdgeCols[1], "white", EdgeCols[2], "white", EdgeCols[3])
                    legend(x = "top", legend = LegendText, lty = c(1,0,1,0,1), col = LegendCols, 
                           bty = "n", cex = LegSize, lwd = LineWidth, horiz = TRUE,
                           inset = LegInset, xpd = NA)
               }
          }
          # Finally, add the axis labels
          mtext(text = xLabel, side = 1, outer = TRUE, line = xLabLine,
                cex = TextSize)
          mtext(text = yLabel, side = 2, outer = TRUE, line = yLabLine,
                cex = TextSize)
     dev.off()
}
rm(FitGenVarData)
