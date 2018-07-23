# This script will create the dispersal evolution graphs
setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/")
load("SimData/DispEvolData.rdata")
SpeedWords <- c("Slow", "Main", "Fast")

# Create objects to hold coordinates for the various arrows to use in the graph
TopArrow <- array(NA, dim = c(3, 2, 2))
TopArrow[1,,] <- c(-75, 600, 16000, 16000)
TopArrow[2,,] <- c(-75, 600, 18000, 18000)
TopArrow[3,,] <- c(-75, 600, 16000, 16000)
SideArrow <- array(NA, dim = c(3, 2, 2))
SideArrow[1,,] <- c(-193, -193, -34000, 10000)
SideArrow[2,,] <- c(-193, -193, -38000, 12000)
SideArrow[3,,] <- c(-193, -193, -34000, 10000)

# Set some graphical parameters for the figures
FigWidth <- 8
FigHeight <- 6
FigMat <- matrix(NA, nrow = 3, ncol = 12)
FigMat[1,] <- c(rep(1,4), rep(2,4), rep(3,4))
FigMat[2,] <- c(rep(4,4), rep(5,4), rep(6,4))
FigMat[3,] <- c(rep(7,4), rep(8,4), rep(9,4))
OuterMar <- c(4, 10, 6, 2)
InnerMar <- c(1.5, 1.5, 1.5, 1.5)
SimSeq <- c(7,8,9,4,5,6,1,2,3)
xRange <- c(-100, 100)
AxisSize <- 1.05
ArrowWidth <- 3
ArrowLength <- 0.25
DispAxisLabel <- expression(paste(Delta, bar(d), sep = ""))
TextSize <- 1.5
xLabLine <- 2
yLabLine <- 1.75
GradLabLine <- 3.5
GradSubLine <- 1.5
AdaptLabLine <- 7
AdaptSubLine <- 5
LowAdj <- 0.05
HighAdj <- 0.95

AllCol <- "lightskyblue"
ExtantCol <- "navyblue"

range(DispEvol[[1]]$DeltaDisp, na.rm = TRUE)
range(DispEvol[[2]]$DeltaDisp)
range(DispEvol[[3]]$DeltaDisp, na.rm = TRUE)
# Set the HistBreaks accordingly...
HistBreaks <- seq(-375, 475, by = 5)

# Now make the figures
for(v in 1:3){
     PlotName <- paste("ResultFigures/", SpeedWords[v], "DispEvol.pdf", sep = "")
     CurData <- DispEvol[[v]]
     ExtantData <- paste("SimData/", SpeedWords[v], "NumExtant.rdata", sep = "")
     load(ExtantData)
     pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
          layout(FigMat)
          par(mar = InnerMar, oma = OuterMar)
          for(i in SimSeq){
               ParamData <- subset(CurData, Param == i)
               ExtantData <- subset(CurData, (Param == i) & (Ext == 1))
               ParamHist <- hist(ParamData$DeltaDisp, plot = FALSE, breaks = HistBreaks)
               plot(ParamHist, col = AllCol, main = "", xlab = "", ylab = "", 
                    las = 1, xlim = xRange, xaxt = "n", yaxt = "n")
               axis(1, cex = AxisSize)
               axis(2, cex.axis = AxisSize)
               hist(ExtantData$DeltaDisp, col = ExtantCol, breaks = HistBreaks,
                    add = TRUE)
               # Add in the number of extant simulations
               SimMessage <- paste("n = ", NumExtant[i], sep = "")
               legend("topleft", legend = SimMessage, lty = 0, text.col = ExtantCol,
                      bty = "n", cex = AxisSize)
               # Add the selection and environmental gradient arrows
               if(i == 7){
                    arrows(x0 = TopArrow[v,1,1], y0 = TopArrow[v,1,2], 
                           x1 = TopArrow[v,2,1], y1 = TopArrow[v,2,2],
                           length = ArrowLength, lwd = ArrowWidth, xpd = NA)
                    arrows(x0 = SideArrow[v,1,1], y0 = SideArrow[v,1,2], 
                           x1 = SideArrow[v,2,1], y1 = SideArrow[v,2,2],
                           length = ArrowLength, lwd = ArrowWidth, xpd = NA)
               }
          }
          # Add the x and y axis labels
          mtext(DispAxisLabel, side = 1, outer = TRUE, line = xLabLine, cex = TextSize)
          mtext("Frequency", side = 2, outer = TRUE, line = yLabLine, cex = TextSize)
          
          # Add the adaptation and gradient text
          mtext("Adaptation potential", side = 2, outer = TRUE, line = AdaptLabLine,
                cex = TextSize)
          mtext("None", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize,
                adj = LowAdj)
          mtext("Low", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize)
          mtext("High", side = 2, outer = TRUE, line = AdaptSubLine, cex = TextSize,
                adj = HighAdj)
          
          mtext("Gradient at range edge", side = 3, outer = TRUE, line = GradLabLine,
                cex = TextSize)
          mtext("Gradual", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
                adj = LowAdj)
          mtext("Moderate", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize)
          mtext("Severe", side = 3, outer = TRUE, line = GradSubLine, cex = TextSize,
                adj = HighAdj)
     dev.off()
}

# Now, test whether the extant simulations are significantly different from
#    a random subset of the total simulations in terms of mean, standard 
#    deviation, and skew 
library(moments)
DispEvolResults <- expand.grid(param = 1:9, speed = 1:3, 
                               value = c("mu", "sigma", "skew"))
DispEvolResults$Obs <- rep(NA, nrow(DispEvolResults))
DispEvolResults$p.value <- rep(NA, nrow(DispEvolResults))
Sims <- 10000
for(v in 1:3){
     for(p in 1:9){
          # Isolate the current data
          AllCurData <- subset(DispEvol[[v]], Param == p)
          ExtantCurData <- subset(AllCurData, Ext == 1)
          # Get the relevant indices from the results data frame
          MuIndex <- which((DispEvolResults$param == p) & (DispEvolResults$speed == v) & 
                                (DispEvolResults$value == "mu"))
          SigmaIndex <- which((DispEvolResults$param == p) & (DispEvolResults$speed == v) & 
                                (DispEvolResults$value == "sigma"))
          SkewIndex <- which((DispEvolResults$param == p) & (DispEvolResults$speed == v) & 
                                (DispEvolResults$value == "skew"))
          # Calculate the observed values
          DispEvolResults$Obs[MuIndex] <- mean(ExtantCurData$DeltaDisp)
          DispEvolResults$Obs[SigmaIndex] <- sd(ExtantCurData$DeltaDisp)
          DispEvolResults$Obs[SkewIndex] <- skewness(ExtantCurData$DeltaDisp)
          # Perform the randomization test for significance
          N <- nrow(ExtantCurData)
          Means <- rep(NA, Sims)
          Sigmas <- rep(NA, Sims)
          Skews <- rep(NA, Sims)
          for(i in 1:Sims){
               ExtantSeq <- sample(1:nrow(AllCurData), size = N, replace = FALSE)
               ExtantSim <- AllCurData[ExtantSeq,]
               Means[i] <- mean(ExtantSim$DeltaDisp)
               Sigmas[i] <- sd(ExtantSim$DeltaDisp)
               Skews[i] <- skewness(ExtantSim$DeltaDisp)
          }
          DispEvolResults$p.value[MuIndex] <- sum(Means >= DispEvolResults$Obs[MuIndex]) / Sims
          DispEvolResults$p.value[SigmaIndex] <- sum(Sigmas <= DispEvolResults$Obs[SigmaIndex]) / Sims
          DispEvolResults$p.value[SkewIndex] <- sum(Skews >= DispEvolResults$Obs[SkewIndex]) / Sims
     }
}

write.csv(DispEvolResults, file = "DispEvolResults.csv", row.names = FALSE, quote = FALSE)
