# This script will create an animated graph of four metrics: 1) mean dispersal 
#    phenotype, 2) mean nich phenotype, 3) genetic variance in dispersal, and 
#    4) genetic variance in the  niche trait. Each metric will be further broken
#    down to highlight extant vs. extinct simulations. The animation will go
#    through time during the range shift in increments of 10 generations.

# Set the working directory and create some graphical objects to be used across
#    all the animations
setwd("~/Desktop/PostdocResearch/ShiftingSlopesOther/")
TextSize <- 1.5
AxisSize <- 1.75
LegSize <- 2
yArrAdjLwr <- 0.2
yArrAdjUpr <- 3
xArrAdj <- c(1.4, 1.45, 1.5)
FigWidth <- 1200
FigHeight <- 900
FigPointSize <- 18
ExtinctLineCol <- rgb(red = 1, green = 0, blue = 0, alpha = 1)
ExtinctShading <- rgb(red = 1, green = 0, blue = 0, alpha = 0.5)
ExtantLineCol <- rgb(red = 0, green = 0, blue = 1, alpha = 1)
ExtantShading <- rgb(red = 0, green = 0, blue = 1, alpha = 0.5)
PrevShading <- gray(level = 0.5, alpha = 0.5)
EdgeSims <- matrix(c(7,8,9,4,5,6,1,2,3), nrow = 3, ncol = 3, byrow = TRUE)
EdgeWords <- c("Gradual edge", "Moderate edge", "Stark edge")
OptWords <- c("Flat", "Shallow", "Steep")
OptLabel <- "Gradient in niche optimum"
xLabel <- "Space"
xLabLine <- 3
yLabLine <- 4
GenSteps <- 10
LengthShift <- 100
GenSeq <- seq(0, LengthShift, by = GenSteps)
LowAdj <- 0.05
HighAdj <- 0.95
OptLabLine <- 10
OptWordLine <- 7.5
OuterMar <- c(5, 12, 4, 2)
InnerMar <- c(0.35,0,0.35,0)
SpeedWords <- c("Slow", "Main", "Fast")
LineWidth <- 2.5

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
          if((nrow(CurData$Extant) > 0) & (nrow(CurData$Extinct) > 0)){
               xLwr[v] <- min(xLwr[v], min(CurData$Extant$x), min(CurData$Extinct$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extant$x), max(CurData$Extinct$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extant$lwr), min(CurData$Extinct$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extant$upr), max(CurData$Extinct$upr))
          } else if(nrow(CurData$Extant) == 0){
               xLwr[v] <- min(xLwr[v], min(CurData$Extinct$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extinct$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extinct$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extinct$upr))
          } else if(nrow(CurData$Extinct) == 0){
               xLwr[v] <- min(xLwr[v], min(CurData$Extant$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extant$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extant$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extant$upr))
          }
     }
}

for(v in 1:3){
     for(i in 1:3){
          CurSims <- EdgeSims[,i]
          # Create a temporary directory to hold the individual images comprising the
          #    animation
          TempDirPath <- "ResultFigures/temp"
          dir.create(TempDirPath)
          # Now loop through each time point to create a still image for that generation
          for(gen in GenSeq){
               if(gen == 0){
                    GenNum <- "000"
               } else if(gen < 100){
                    GenNum <- paste("0", gen, sep = "")
               } else{
                    GenNum <- gen
               }
               FigName <- paste(TempDirPath, "/gen", GenNum, ".png", sep = "")
               png(filename = FigName, width = FigWidth, height = FigHeight, pointsize = FigPointSize)
                    par(mfrow = c(3,1), oma = OuterMar, mar = InnerMar)
                    for(j in CurSims){
                         CurExtant <- subset(DispData[[j]][[v]]$Extant, g == gen)
                         CurExtinct <- subset(DispData[[j]][[v]]$Extinct, g == gen)
                         CurExtant <- CurExtant[order(CurExtant$x),]
                         CurExtinct <- CurExtinct[order(CurExtinct$x),]
                         if(j == CurSims[3]){
                              plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                                   main = "", xlab = "", ylab = "", cex.axis = AxisSize, las = 1)
                         }else{
                              plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                                   main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize,
                                   las = 1)
                              axis(1, labels = FALSE)
                         }
                         # Determine the sequence of all previous generations and add them
                         #    to the plot
                         CurIndex <- which(GenSeq == gen)
                         PrevGens <- GenSeq[1:CurIndex]
                         if(length(PrevGens) > 1){
                              for(k in PrevGens){
                                   PreviousExtant <- subset(DispData[[j]][[v]]$Extant, g == k)
                                   PreviousExtinct <- subset(DispData[[j]][[v]]$Extinct, g == k)
                                   PreviousExtant <- PreviousExtant[order(PreviousExtant$x),]
                                   PreviousExtinct <- PreviousExtinct[order(PreviousExtinct$x),]
                                   if(nrow(PreviousExtant) > 0){
                                        lines(x = PreviousExtant$x, y = PreviousExtant$dBar, lty = 1,
                                              col = "black")
                                        polygon(x = c(PreviousExtant$x, PreviousExtant$x[nrow(PreviousExtant):1]),
                                                y = c(PreviousExtant$lwr, PreviousExtant$upr[nrow(PreviousExtant):1]),
                                                border = NA, col = PrevShading)
                                   }
                                   if(nrow(PreviousExtinct) > 0){
                                        lines(x = PreviousExtinct$x, y = PreviousExtinct$dBar, lty = 2,
                                              col = "black")
                                        polygon(x = c(PreviousExtinct$x, PreviousExtinct$x[nrow(PreviousExtinct):1]),
                                                y = c(PreviousExtinct$lwr, PreviousExtinct$upr[nrow(PreviousExtinct):1]),
                                                border = NA, col = PrevShading)
                                   }
                              }
                         }
                         # Now add the current generation to the plot
                         if(nrow(CurExtant) > 0){
                              lines(x = CurExtant$x, y = CurExtant$dBar, lty = 1,
                                    col = ExtantLineCol, lwd = LineWidth)
                              polygon(x = c(CurExtant$x, CurExtant$x[nrow(CurExtant):1]),
                                      y = c(CurExtant$lwr, CurExtant$upr[nrow(CurExtant):1]),
                                      border = NA, col = ExtantShading)
                         }
                         if(nrow(CurExtinct) > 0){
                              lines(x = CurExtinct$x, y = CurExtinct$dBar, lty = 2,
                                    col = ExtinctLineCol, lwd = LineWidth)
                              polygon(x = c(CurExtinct$x, CurExtinct$x[nrow(CurExtinct):1]),
                                      y = c(CurExtinct$lwr, CurExtinct$upr[nrow(CurExtinct):1]),
                                      border = NA, col = ExtinctShading)
                         }
                    }
                    # Add a legend on the last plot
                    legend(x = "topleft", legend = c("Extant", "Extinct"), lty = c(1,2),
                           col = c(ExtantLineCol, ExtinctLineCol), bty = "n", cex = LegSize,
                           lwd = LineWidth)
                    # Add the edge description and generation
                    mtext(text = EdgeWords[i], side = 3, outer = TRUE, cex = TextSize)
                    mtext(text = paste("Generation", gen, sep = " "), side = 3,
                          outer = TRUE, cex = TextSize, adj = 0.95)
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
          # Create the gif in the main directory and then remove the temporary directory
          GifName <- paste("ResultFigures/", SpeedWords[v], strsplit(x = EdgeWords[i], split = " ")[[1]][1], "MuDisp.gif", sep = "")
          SysCommand1 <- paste("convert -delay 50 ", TempDirPath, "/gen*.png ", GifName, sep = "")
          SysCommand2 <- paste("rm -r ", TempDirPath, sep = "")
          system(SysCommand1)
          system(SysCommand2)
     }
}
rm(DispData)

# ------------------------------------------------------ 2) mean niche phenotype
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
          if((nrow(CurData$Extant) > 0) & (nrow(CurData$Extinct) > 0)){
               xLwr[v] <- min(xLwr[v], min(CurData$Extant$x), min(CurData$Extinct$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extant$x), max(CurData$Extinct$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extant$lwr), min(CurData$Extinct$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extant$upr), max(CurData$Extinct$upr))
          } else if(nrow(CurData$Extant) == 0){
               xLwr[v] <- min(xLwr[v], min(CurData$Extinct$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extinct$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extinct$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extinct$upr))
          } else if(nrow(CurData$Extinct) == 0){
               xLwr[v] <- min(xLwr[v], min(CurData$Extant$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extant$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extant$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extant$upr))
          }
     }
}

for(v in 1:3){
     for(i in 1:3){
          CurSims <- EdgeSims[,i]
          # Create a temporary directory to hold the individual images comprising the
          #    animation
          TempDirPath <- "ResultFigures/temp"
          dir.create(TempDirPath)
          # Now loop through each time point to create a still image for that generation
          for(gen in GenSeq){
               if(gen == 0){
                    GenNum <- "000"
               } else if(gen < 100){
                    GenNum <- paste("0", gen, sep = "")
               } else{
                    GenNum <- gen
               }
               FigName <- paste(TempDirPath, "/gen", GenNum, ".png", sep = "")
               png(filename = FigName, width = FigWidth, height = FigHeight, pointsize = FigPointSize)
                    par(mfrow = c(3,1), oma = OuterMar, mar = InnerMar)
                    for(j in CurSims){
                         CurExtant <- subset(FitData[[j]][[v]]$Extant, g == gen)
                         CurExtinct <- subset(FitData[[j]][[v]]$Extinct, g == gen)
                         CurExtant <- CurExtant[order(CurExtant$x),]
                         CurExtinct <- CurExtinct[order(CurExtinct$x),]
                         if(j == CurSims[3]){
                              plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                                   main = "", xlab = "", ylab = "", cex.axis = AxisSize, las = 1)
                         }else{
                              plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                                   main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize,
                                   las = 1)
                              axis(1, labels = FALSE)
                         }
                         # Determine the sequence of all previous generations and add them
                         #    to the plot
                         CurIndex <- which(GenSeq == gen)
                         PrevGens <- GenSeq[1:CurIndex]
                         if(length(PrevGens) > 1){
                              for(k in PrevGens){
                                   PreviousExtant <- subset(FitData[[j]][[v]]$Extant, g == k)
                                   PreviousExtinct <- subset(FitData[[j]][[v]]$Extinct, g == k)
                                   PreviousExtant <- PreviousExtant[order(PreviousExtant$x),]
                                   PreviousExtinct <- PreviousExtinct[order(PreviousExtinct$x),]
                                   if(nrow(PreviousExtant) > 0){
                                        lines(x = PreviousExtant$x, y = PreviousExtant$wBar, lty = 1,
                                              col = "black")
                                        polygon(x = c(PreviousExtant$x, PreviousExtant$x[nrow(PreviousExtant):1]),
                                                y = c(PreviousExtant$lwr, PreviousExtant$upr[nrow(PreviousExtant):1]),
                                                border = NA, col = PrevShading)
                                   }
                                   if(nrow(PreviousExtinct) > 0){
                                        lines(x = PreviousExtinct$x, y = PreviousExtinct$wBar, lty = 2,
                                              col = "black")
                                        polygon(x = c(PreviousExtinct$x, PreviousExtinct$x[nrow(PreviousExtinct):1]),
                                                y = c(PreviousExtinct$lwr, PreviousExtinct$upr[nrow(PreviousExtinct):1]),
                                                border = NA, col = PrevShading)
                                   }
                              }
                         }
                         # Now add the current generation to the plot
                         if(nrow(CurExtant) > 0){
                              lines(x = CurExtant$x, y = CurExtant$wBar, lty = 1,
                                    col = ExtantLineCol, lwd = LineWidth)
                              polygon(x = c(CurExtant$x, CurExtant$x[nrow(CurExtant):1]),
                                      y = c(CurExtant$lwr, CurExtant$upr[nrow(CurExtant):1]),
                                      border = NA, col = ExtantShading)
                         }
                         if(nrow(CurExtinct) > 0){
                              lines(x = CurExtinct$x, y = CurExtinct$wBar, lty = 2,
                                    col = ExtinctLineCol, lwd = LineWidth)
                              polygon(x = c(CurExtinct$x, CurExtinct$x[nrow(CurExtinct):1]),
                                      y = c(CurExtinct$lwr, CurExtinct$upr[nrow(CurExtinct):1]),
                                      border = NA, col = ExtinctShading)
                         }
                    }
                    # Add a legend on the last plot
                    legend(x = "topleft", legend = c("Extant", "Extinct"), lty = c(1,2),
                           col = c(ExtantLineCol, ExtinctLineCol), bty = "n", cex = LegSize,
                           lwd = LineWidth)
                    # Add the edge description and generation
                    mtext(text = EdgeWords[i], side = 3, outer = TRUE, cex = TextSize)
                    mtext(text = paste("Generation", gen, sep = " "), side = 3,
                          outer = TRUE, cex = TextSize, adj = 0.95)
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
          # Create the gif in the main directory and then remove the temporary directory
          GifName <- paste("ResultFigures/", SpeedWords[v], strsplit(x = EdgeWords[i], split = " ")[[1]][1], "MuFit.gif", sep = "")
          SysCommand1 <- paste("convert -delay 50 ", TempDirPath, "/gen*.png ", GifName, sep = "")
          SysCommand2 <- paste("rm -r ", TempDirPath, sep = "")
          system(SysCommand1)
          system(SysCommand2)
     }
}
rm(FitData)

# --------------------------------------------- 3) genetic variance in dispersal
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
          if((nrow(CurData$Extant) > 0) & (nrow(CurData$Extinct) > 0)){
               xLwr[v] <- min(xLwr[v], min(CurData$Extant$x), min(CurData$Extinct$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extant$x), max(CurData$Extinct$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extant$lwr), min(CurData$Extinct$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extant$upr), max(CurData$Extinct$upr))
          } else if(nrow(CurData$Extant) == 0){
               xLwr[v] <- min(xLwr[v], min(CurData$Extinct$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extinct$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extinct$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extinct$upr))
          } else if(nrow(CurData$Extinct) == 0){
               xLwr[v] <- min(xLwr[v], min(CurData$Extant$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extant$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extant$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extant$upr))
          }
     }
}

for(v in 1:3){
     for(i in 1:3){
          CurSims <- EdgeSims[,i]
          # Create a temporary directory to hold the individual images comprising the
          #    animation
          TempDirPath <- "ResultFigures/temp"
          dir.create(TempDirPath)
          # Now loop through each time point to create a still image for that generation
          for(gen in GenSeq){
               if(gen == 0){
                    GenNum <- "000"
               } else if(gen < 100){
                    GenNum <- paste("0", gen, sep = "")
               } else{
                    GenNum <- gen
               }
               FigName <- paste(TempDirPath, "/gen", GenNum, ".png", sep = "")
               png(filename = FigName, width = FigWidth, height = FigHeight, pointsize = FigPointSize)
                    par(mfrow = c(3,1), oma = OuterMar, mar = InnerMar)
                    for(j in CurSims){
                         CurExtant <- subset(DispGenVarData[[j]][[v]]$Extant, g == gen)
                         CurExtinct <- subset(DispGenVarData[[j]][[v]]$Extinct, g == gen)
                         CurExtant <- CurExtant[order(CurExtant$x),]
                         CurExtinct <- CurExtinct[order(CurExtinct$x),]
                         if(j == CurSims[3]){
                              plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                                   main = "", xlab = "", ylab = "", cex.axis = AxisSize, las = 1)
                         }else{
                              plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                                   main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize,
                                   las = 1)
                              axis(1, labels = FALSE)
                         }
                         # Determine the sequence of all previous generations and add them
                         #    to the plot
                         CurIndex <- which(GenSeq == gen)
                         PrevGens <- GenSeq[1:CurIndex]
                         if(length(PrevGens) > 1){
                              for(k in PrevGens){
                                   PreviousExtant <- subset(DispGenVarData[[j]][[v]]$Extant, g == k)
                                   PreviousExtinct <- subset(DispGenVarData[[j]][[v]]$Extinct, g == k)
                                   PreviousExtant <- PreviousExtant[order(PreviousExtant$x),]
                                   PreviousExtinct <- PreviousExtinct[order(PreviousExtinct$x),]
                                   if(nrow(PreviousExtant) > 0){
                                        lines(x = PreviousExtant$x, y = PreviousExtant$GenVar, lty = 1,
                                              col = "black")
                                        polygon(x = c(PreviousExtant$x, PreviousExtant$x[nrow(PreviousExtant):1]),
                                                y = c(PreviousExtant$lwr, PreviousExtant$upr[nrow(PreviousExtant):1]),
                                                border = NA, col = PrevShading)
                                   }
                                   if(nrow(PreviousExtinct) > 0){
                                        lines(x = PreviousExtinct$x, y = PreviousExtinct$GenVar, lty = 2,
                                              col = "black")
                                        polygon(x = c(PreviousExtinct$x, PreviousExtinct$x[nrow(PreviousExtinct):1]),
                                                y = c(PreviousExtinct$lwr, PreviousExtinct$upr[nrow(PreviousExtinct):1]),
                                                border = NA, col = PrevShading)
                                   }
                              }
                         }
                         # Now add the current generation to the plot
                         if(nrow(CurExtant) > 0){
                              lines(x = CurExtant$x, y = CurExtant$GenVar, lty = 1,
                                    col = ExtantLineCol, lwd = LineWidth)
                              polygon(x = c(CurExtant$x, CurExtant$x[nrow(CurExtant):1]),
                                      y = c(CurExtant$lwr, CurExtant$upr[nrow(CurExtant):1]),
                                      border = NA, col = ExtantShading)
                         }
                         if(nrow(CurExtinct) > 0){
                              lines(x = CurExtinct$x, y = CurExtinct$GenVar, lty = 2,
                                    col = ExtinctLineCol, lwd = LineWidth)
                              polygon(x = c(CurExtinct$x, CurExtinct$x[nrow(CurExtinct):1]),
                                      y = c(CurExtinct$lwr, CurExtinct$upr[nrow(CurExtinct):1]),
                                      border = NA, col = ExtinctShading)
                         }
                    }
                    # Add a legend on the last plot
                    legend(x = "topleft", legend = c("Extant", "Extinct"), lty = c(1,2),
                           col = c(ExtantLineCol, ExtinctLineCol), bty = "n", cex = LegSize,
                           lwd = LineWidth)
                    # Add the edge description and generation
                    mtext(text = EdgeWords[i], side = 3, outer = TRUE, cex = TextSize)
                    mtext(text = paste("Generation", gen, sep = " "), side = 3,
                          outer = TRUE, cex = TextSize, adj = 0.95)
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
                          cex = TextSize)
               dev.off()
          }
          # Create the gif in the main directory and then remove the temporary directory
          GifName <- paste("ResultFigures/", SpeedWords[v], strsplit(x = EdgeWords[i], split = " ")[[1]][1], "DispGenVar.gif", sep = "")
          SysCommand1 <- paste("convert -delay 50 ", TempDirPath, "/gen*.png ", GifName, sep = "")
          SysCommand2 <- paste("rm -r ", TempDirPath, sep = "")
          system(SysCommand1)
          system(SysCommand2)
     }
}
rm(DispGenVarData)

# --------------------------------------------- 3) genetic variance in dispersal
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
          if((nrow(CurData$Extant) > 0) & (nrow(CurData$Extinct) > 0)){
               xLwr[v] <- min(xLwr[v], min(CurData$Extant$x), min(CurData$Extinct$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extant$x), max(CurData$Extinct$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extant$lwr), min(CurData$Extinct$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extant$upr), max(CurData$Extinct$upr))
          } else if(nrow(CurData$Extant) == 0){
               xLwr[v] <- min(xLwr[v], min(CurData$Extinct$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extinct$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extinct$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extinct$upr))
          } else if(nrow(CurData$Extinct) == 0){
               xLwr[v] <- min(xLwr[v], min(CurData$Extant$x))
               xUpr[v] <- max(xUpr[v], max(CurData$Extant$x))
               yLwr[v] <- min(yLwr[v], min(CurData$Extant$lwr))
               yUpr[v] <- max(yUpr[v], max(CurData$Extant$upr))
          }
     }
}

for(v in 1:3){
     for(i in 1:3){
          CurSims <- EdgeSims[,i]
          # Create a temporary directory to hold the individual images comprising the
          #    animation
          TempDirPath <- "ResultFigures/temp"
          dir.create(TempDirPath)
          # Now loop through each time point to create a still image for that generation
          for(gen in GenSeq){
               if(gen == 0){
                    GenNum <- "000"
               } else if(gen < 100){
                    GenNum <- paste("0", gen, sep = "")
               } else{
                    GenNum <- gen
               }
               FigName <- paste(TempDirPath, "/gen", GenNum, ".png", sep = "")
               png(filename = FigName, width = FigWidth, height = FigHeight, pointsize = FigPointSize)
               par(mfrow = c(3,1), oma = OuterMar, mar = InnerMar)
               for(j in CurSims){
                    CurExtant <- subset(FitGenVarData[[j]][[v]]$Extant, g == gen)
                    CurExtinct <- subset(FitGenVarData[[j]][[v]]$Extinct, g == gen)
                    CurExtant <- CurExtant[order(CurExtant$x),]
                    CurExtinct <- CurExtinct[order(CurExtinct$x),]
                    if(j == CurSims[3]){
                         plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                              main = "", xlab = "", ylab = "", cex.axis = AxisSize, las = 1)
                    }else{
                         plot(x = NA, y = NA, xlim = c(xLwr[v], xUpr[v]), ylim = c(yLwr[v], yUpr[v]),
                              main = "", xlab = "", ylab = "", xaxt = "n", cex.axis = AxisSize,
                              las = 1)
                         axis(1, labels = FALSE)
                    }
                    # Determine the sequence of all previous generations and add them
                    #    to the plot
                    CurIndex <- which(GenSeq == gen)
                    PrevGens <- GenSeq[1:CurIndex]
                    if(length(PrevGens) > 1){
                         for(k in PrevGens){
                              PreviousExtant <- subset(FitGenVarData[[j]][[v]]$Extant, g == k)
                              PreviousExtinct <- subset(FitGenVarData[[j]][[v]]$Extinct, g == k)
                              PreviousExtant <- PreviousExtant[order(PreviousExtant$x),]
                              PreviousExtinct <- PreviousExtinct[order(PreviousExtinct$x),]
                              if(nrow(PreviousExtant) > 0){
                                   lines(x = PreviousExtant$x, y = PreviousExtant$GenVar, lty = 1,
                                         col = "black")
                                   polygon(x = c(PreviousExtant$x, PreviousExtant$x[nrow(PreviousExtant):1]),
                                           y = c(PreviousExtant$lwr, PreviousExtant$upr[nrow(PreviousExtant):1]),
                                           border = NA, col = PrevShading)
                              }
                              if(nrow(PreviousExtinct) > 0){
                                   lines(x = PreviousExtinct$x, y = PreviousExtinct$GenVar, lty = 2,
                                         col = "black")
                                   polygon(x = c(PreviousExtinct$x, PreviousExtinct$x[nrow(PreviousExtinct):1]),
                                           y = c(PreviousExtinct$lwr, PreviousExtinct$upr[nrow(PreviousExtinct):1]),
                                           border = NA, col = PrevShading)
                              }
                         }
                    }
                    # Now add the current generation to the plot
                    if(nrow(CurExtant) > 0){
                         lines(x = CurExtant$x, y = CurExtant$GenVar, lty = 1,
                               col = ExtantLineCol, lwd = LineWidth)
                         polygon(x = c(CurExtant$x, CurExtant$x[nrow(CurExtant):1]),
                                 y = c(CurExtant$lwr, CurExtant$upr[nrow(CurExtant):1]),
                                 border = NA, col = ExtantShading)
                    }
                    if(nrow(CurExtinct) > 0){
                         lines(x = CurExtinct$x, y = CurExtinct$GenVar, lty = 2,
                               col = ExtinctLineCol, lwd = LineWidth)
                         polygon(x = c(CurExtinct$x, CurExtinct$x[nrow(CurExtinct):1]),
                                 y = c(CurExtinct$lwr, CurExtinct$upr[nrow(CurExtinct):1]),
                                 border = NA, col = ExtinctShading)
                    }
               }
               # Add a legend on the last plot
               legend(x = "topleft", legend = c("Extant", "Extinct"), lty = c(1,2),
                      col = c(ExtantLineCol, ExtinctLineCol), bty = "n", cex = LegSize,
                      lwd = LineWidth)
               # Add the edge description and generation
               mtext(text = EdgeWords[i], side = 3, outer = TRUE, cex = TextSize)
               mtext(text = paste("Generation", gen, sep = " "), side = 3,
                     outer = TRUE, cex = TextSize, adj = 0.95)
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
                     cex = TextSize)
               dev.off()
          }
          # Create the gif in the main directory and then remove the temporary directory
          GifName <- paste("ResultFigures/", SpeedWords[v], strsplit(x = EdgeWords[i], split = " ")[[1]][1], "FitGenVar.gif", sep = "")
          SysCommand1 <- paste("convert -delay 50 ", TempDirPath, "/gen*.png ", GifName, sep = "")
          SysCommand2 <- paste("rm -r ", TempDirPath, sep = "")
          system(SysCommand1)
          system(SysCommand2)
     }
}
rm(FitGenVarData)
