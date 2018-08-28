# Make a figure displaying the range shift of an individual simulation

setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/SchematicFigures/")
source("6278_Sim1/parameters.R")
source("~/Desktop/RangeShifts/ShiftingSlopesCode/SimFunctions.R")
SumStats <- read.csv("6278_Sim1/SummaryStats.csv")

# Define the time points that will be illustrated and extract just those data
TimePoints <- seq(50, 90, by = 5)
FigData <- subset(SumStats, (gen %in% TimePoints), select = c(gen, beta, x, y, 
                                                              abund, muFit, muDisp))

# Create the large array holding the matrices that will hold the population 
#    values at different time points
NumCols <- width * length(TimePoints)
NumRows <- max(FigData$x) - min(FigData$x)
FigMats <- array(NA, dim = c(3, NumRows, NumCols))

# Now populate the figure matrices
xCorrection <- abs(min(FigData$x))
yCorrections <- ((1:length(TimePoints)) - 1) * width
Betas <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, v = 1, eta = eta)
for(t in 1:length(TimePoints)){
     CurData <- subset(FigData, gen == TimePoints[t])
     if(t == 1){
          CurBeta <- BetaInit
     } else{
          CurBeta <- Betas[TimePoints[t] - 50]
     }
     for(i in 1:dim(CurData)[1]){
          FigMats[1,CurData$x[i] + xCorrection,CurData$y[i] + yCorrections[t]] <- 
               CurData$abund[i]
          FigMats[2,CurData$x[i] + xCorrection,CurData$y[i] + yCorrections[t]] <- 
               log(CurData$muDisp[i], base = 10)
          Zopt <- lambda * (CurData$x[i] * eta - CurBeta)
          RelFit <- exp(-1 * (CurData$muFit[i] - Zopt)^2 / (2*omega^2))
          FigMats[3,CurData$x[i] + xCorrection,CurData$y[i] + yCorrections[t]] <- 
               RelFit
     }
}

# Set some graphical parameters to use for the subsequent figures
InnerMar <- c(2, 4, 3, 2) + 0.1
ColFunc <- colorRampPalette(c("lightyellow", "darkorange", "darkred"))
TextSize <- 1.15
FigWidth <- 8
FigHeight <- 2.5
library(plot3D)
LineSeq <- ((1:length(TimePoints)) * width) / dim(FigMats[1,,])[2]
LineWidth <- 0.5
TimeLabels <- paste("t = ", TimePoints-50, sep = "")
TimeSize <- 1
PlotName <- "SimExample.pdf"
pdf(file = PlotName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     par(mfrow = c(1,3), mar = InnerMar)
     # Abundance
     image2D(z = FigMats[1,,], xlab = "", ylab = "", main = "", xaxt = "n", 
             yaxt = "n", NAcol = "black", col = ColFunc(10000))
     mtext("Abundance", side = 3, line = 0.5, cex = TextSize)
     abline(h = LineSeq[1:8], lty = 1, lwd = LineWidth, col = "white")
     text(x = -0.125, y = LineSeq - 0.05, labels = TimeLabels, cex = TimeSize, xpd = NA)
     text(x = 0.05, y = LineSeq[9] - 0.05, labels = expression(bold("a")), cex = TextSize, col = "white")
     # Dispersal
     image2D(z = FigMats[2,,], xlab = "", ylab = "", main = "", xaxt = "n", 
             yaxt = "n", NAcol = "black", col = ColFunc(10000))
     mtext("Dispersal", side = 3, line = 0.5, cex = TextSize)
     abline(h = LineSeq[1:8], lty = 1, lwd = LineWidth, col = "white")
     text(x = -0.125, y = LineSeq - 0.05, labels = TimeLabels, cex = TimeSize, xpd = NA)
     text(x = 0.05, y = LineSeq[9] - 0.05, labels = expression(bold("b")), cex = TextSize, col = "white")
     # Fitness
     image2D(z = FigMats[3,,], xlab = "", ylab = "", main = "", xaxt = "n", 
             yaxt = "n", NAcol = "black", col = ColFunc(10000))
     mtext("Fitness", side = 3, line = 0.5, cex = TextSize)
     abline(h = LineSeq[1:8], lty = 1, lwd = LineWidth, col = "white")
     text(x = -0.125, y = LineSeq - 0.05, labels = TimeLabels, cex = TimeSize, xpd = NA)
     text(x = 0.05, y = LineSeq[9] - 0.05, labels = expression(bold("c")), cex = TextSize, col = "white")
dev.off()
