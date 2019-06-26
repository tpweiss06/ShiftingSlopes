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
LegInset <- -0.6
LetterSize <- 2.5
OuterMar <- c(11, 9, 2, 2)

# -------------------------------------------------- 1) mean dispersal phenotype
yLabel <- expression(bar(d))
load("SimData/DispDataNew.rdata")
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
               # Add a legend
               if(n == 3){
                    legend(x = "bottom", legend = EdgeWords, lty = 1, col = EdgeCols, 
                           cex = LegSize, lwd = LineWidth, horiz = TRUE,
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
yLabel <- "Genetic variance in dispersal"
load("SimData/DispGenVarDataNew.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- rep(Inf, 3)
xUpr <- rep(-Inf, 3)
yLwr <- rep(Inf, 3)
yUpr <- rep(-Inf, 3)
for(p in 1:9){
     for(v in 1:3){
          CurData <- DispGenVarData[[p]][[v]]
          xLwr[v] <- min(xLwr[v], min(CurData$x))
          xUpr[v] <- max(xUpr[v], max(CurData$x))
          yLwr[v] <- min(yLwr[v], min(CurData$lwr))
          yUpr[v] <- max(yUpr[v], max(CurData$upr))
     }
}

for(v in 1:3){
     FigName <- paste("ResultFigures/", SpeedWords[v], "DispGenVar.pdf", sep = "")
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
               # Add a legend
               if(n == 3){
                    legend(x = "bottom", legend = EdgeWords, lty = 1, col = EdgeCols, 
                           cex = LegSize, lwd = LineWidth, horiz = TRUE,
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
LegInset <- -0.525
LetterSize <- 2
OuterMar <- c(9, 6, 2, 2)
#############

# ----------------------------------------------------- 3) mean relative fitness
yLabel <- expression(bar(w))
load("SimData/FitDataNew.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- rep(Inf, 3)
xUpr <- rep(-Inf, 3)
yLwr <- rep(Inf, 3)
yUpr <- rep(-Inf, 3)
for(p in 1:9){
     for(v in 1:3){
          CurData <- FitData[[p]][[v]]
          xLwr[v] <- min(xLwr[v], min(CurData$x))
          xUpr[v] <- max(xUpr[v], max(CurData$x))
          yLwr[v] <- min(yLwr[v], min(CurData$lwr))
          yUpr[v] <- max(yUpr[v], max(CurData$upr))
     }
}

for(v in 1:3){
     FigName <- paste("ResultFigures/", SpeedWords[v], "MeanFit.pdf", sep = "")
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
                         CurData <- subset(FitData[[CurParam]][[v]], g == PlotGens[i])
                         CurData <- CurData[order(CurData$x),]
                         lines(x = CurData$x, y = CurData$wBar, lty = 1,
                               col = EdgeCols[e], lwd = LineWidth)
                         lines(x = CurData$x, y = CurData$lwr, lty = 3,
                               col = EdgeCols[e])
                         lines(x = CurData$x, y = CurData$upr, lty = 3,
                               col = EdgeCols[e])
                         if(e == 1){
                              text(x = BetaPos[v,i], y = (1-GenTextAdj)*yUpr[v], 
                                   labels = paste("t = ", PlotGens[i], sep = ""),
                                   cex = LetterSize)
                         }
                    }
               }
               # Add a letter to the current plot
               text(x = xLetterAdj * xLwr[v], y = yLetterAdj * yUpr[v], 
                    labels = LetterSeq[n], cex = LetterSize)
               # Add a legend
               if(n == 2){
                    legend(x = "bottom", legend = EdgeWords, lty = 1, col = EdgeCols, 
                           cex = LegSize, lwd = LineWidth, horiz = TRUE,
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
rm(FitData)

#-------------------------------------------------- 4) Graph the niche genotypes
yLabel <- "Niche genotypes"
load("SimData/MismatchData.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- rep(Inf, 3)
xUpr <- rep(-Inf, 3)
yLwr <- rep(Inf, 3)
yUpr <- rep(-Inf, 3)
for(p in 1:9){
     for(v in 1:3){
          CurData <- MismatchData[[p]][[v]]
          xLwr[v] <- min(xLwr[v], min(CurData$x))
          xUpr[v] <- max(xUpr[v], max(CurData$x))
          yLwr[v] <- min(yLwr[v], min(CurData$lwr))
          yUpr[v] <- max(yUpr[v], max(CurData$upr))
     }
}

for(v in 1:3){
     FigName <- paste("ResultFigures/", SpeedWords[v], "MeanNiche.pdf", sep = "")
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
               # Add a legend
               if(n == 2){
                    legend(x = "bottom", legend = EdgeWords, lty = 1, col = EdgeCols, 
                           cex = LegSize, lwd = LineWidth, horiz = TRUE,
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

# ----------------------------------------------- 5) genetic variance in fitness
yLabel <- "Genetic variance in fitness"
load("SimData/FitGenVarDataNew.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- rep(Inf, 3)
xUpr <- rep(-Inf, 3)
yLwr <- rep(Inf, 3)
yUpr <- rep(-Inf, 3)
for(p in 1:9){
     for(v in 1:3){
          CurData <- FitGenVarData[[p]][[v]]
          xLwr[v] <- min(xLwr[v], min(CurData$x))
          xUpr[v] <- max(xUpr[v], max(CurData$x))
          yLwr[v] <- min(yLwr[v], min(CurData$lwr))
          yUpr[v] <- max(yUpr[v], max(CurData$upr))
     }
}

for(v in 1:3){
     FigName <- paste("ResultFigures/", SpeedWords[v], "FitGenVar.pdf", sep = "")
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
               # Add a legend
               if(n == 2){
                    legend(x = "bottom", legend = EdgeWords, lty = 1, col = EdgeCols, 
                           cex = LegSize, lwd = LineWidth, horiz = TRUE,
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

############# Now update some graphical parameters for the extinction graph
yLabLine <- 4
TextSize <- 1.5
AxisSize <- 2
LegSize <- 2.25
LegInset <- -0.55
LetterSize <- 2
LetterPos <- c(2, 0.925)
LetterSeq <- c("c", "b", "a")
OuterMar <- c(6,2,0,0)
FigWidth <- 12
FigHeight <- 4.5
xRange <- c(0,100)
yRange <- c(0,1)
TimeAxisSeq1 <- seq(0, 100, by = 20)
TimeAxisSea2 <- seq(0, 100, by = 5)
ExtAxisSeq1 <- seq(0, 1, by = 0.2)
ExtAxisSeq2 <- seq(0, 1, by = 0.05)
ExtAxisLabel <- "Extinction Probability"
#############

# ---------------------------------------------------- 6) Extinction probability
# This section will remake the extinction graph with the same color scheme
#    as the new graphs
# Load the extinction data
load("SimData/Extinctions.rdata")
NumSims <- 200
for(i in 1:3){
     for(p in 1:9){
          Extinctions[i,p,] <- cumsum(Extinctions[i,p,]) / NumSims
     }
}

# Make the figures
for(i in 1:3){
     PlotName <- paste("ResultFigures/", SpeedWords[i], "Extinction.pdf", sep = "")
     pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          par(mfrow = c(1,3), oma = OuterMar)
          for(n in 3:1){
               plot(x = NA, y = NA, xlim = xRange, ylim = yRange, xlab = "",
                    ylab = "", main = "", las = 1, xaxt = "n", yaxt = "n")
               for(j in 1:3){
                    lines(x = 1:100, y = Extinctions[i,SimMat[n,j],], lwd = LineWidth,
                          col = EdgeCols[j])
               }
               # Add the letter and the horizontal line at 1
               text(x = LetterPos[1], y = LetterPos[2], labels = LetterSeq[n], 
                    cex = LetterSize)
               abline(h = 1, lty = 1, col = "grey")
               # Add the axes
               axis(1, at = TimeAxisSeq1, cex.axis = AxisSize)
               axis(1, at = TimeAxisSea2, labels = FALSE, tcl = -0.25)
               axis(2, at = ExtAxisSeq1, las = 1, cex.axis = AxisSize)
               axis(2, at = ExtAxisSeq2, labels = FALSE, tcl = -0.25)
               # Add a legend
               if(n == 2){
                    legend(x = "bottom", legend = EdgeWords, lty = 1, col = EdgeCols, 
                           cex = LegSize, lwd = LineWidth, horiz = TRUE,
                           inset = LegInset, xpd = NA)
               }
          }
          # Add the x and y axis labels
          mtext("Generation", side = 1, outer = TRUE, line = -1, cex = TextSize)
          mtext(ExtAxisLabel, side = 2, outer = TRUE, line = -0.25, cex = TextSize)
     dev.off()
}


