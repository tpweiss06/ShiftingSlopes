# This script will create an animated graph of four metrics: 1) mean dispersal 
#    phenotype, 2) mean nich phenotype, 3) genetic variance in dispersal, and 
#    4) genetic variance in the  niche trait. The animation will go
#    through time during the range shift in increments of 10 generations.

# Set the working directory and create some graphical objects to be used across
#    all the animations
setwd("~/Desktop/PostdocResearch/ShiftingSlopesOther/")
FigWidth <- 1200
FigHeight <- 900
FigPointSize <- 18
GenSteps <- 10
LengthShift <- 2000
GenSeq <- seq(1510, LengthShift, by = GenSteps)
OuterMar <- c(2, 2, 2, 2)
InnerMar <- c(0,0,0,0)
LineWidth <- 2.5

# -------------------------------------------------- 1) mean dispersal phenotype
load("SimData/EquilibriumDispData.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- Inf
xUpr <- -Inf
yLwr <- Inf
yUpr <- -Inf
for(p in 1:9){
          CurData <- DispData[[p]]
          xLwr <- min(xLwr, min(CurData$x))
          xUpr <- max(xUpr, max(CurData$x))
          yLwr <- min(yLwr, min(CurData$lwr))
          yUpr <- max(yUpr, max(CurData$upr))
}

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
          par(mfrow = c(3,3), oma = OuterMar, mar = InnerMar)
          for(i in 1:9){
               CurData <- subset(DispData[[i]], g == gen)
               CurData <- CurData[order(CurData$x),]
               plot(x = NA, y = NA, xlim = c(xLwr, xUpr), ylim = c(yLwr, yUpr),
                    main = "", xlab = "", ylab = "", axes = FALSE)
               lines(x = CurData$x, y = CurData$dBar, lwd = LineWidth)
               lines(x = CurData$x, y = CurData$lwr, lty = 3)
               lines(x = CurData$x, y = CurData$upr, lty = 3)
               box()
          }
     dev.off()
}
# Create the gif in the main directory and then remove the temporary directory
GifName <- "ResultFigures/Animations/Equilibrium/MuDisp.gif"
SysCommand1 <- paste("convert -delay 50 ", TempDirPath, "/gen*.png ", GifName, sep = "")
SysCommand2 <- paste("rm -r ", TempDirPath, sep = "")
system(SysCommand1)
system(SysCommand2)
rm(DispData)

# ------------------------------------------------------ 2) mean niche phenotype
load("SimData/EquilibriumMismatchData.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- Inf
xUpr <- -Inf
yLwr <- Inf
yUpr <- -Inf
for(p in 1:9){
     CurData <- MismatchData[[p]]
     xLwr <- min(xLwr, min(CurData$x))
     xUpr <- max(xUpr, max(CurData$x))
     yLwr <- min(yLwr, min(CurData$lwr))
     yUpr <- max(yUpr, max(CurData$upr))
}

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
     par(mfrow = c(3,3), oma = OuterMar, mar = InnerMar)
     for(i in 1:9){
          CurData <- subset(MismatchData[[i]], g == gen)
          CurData <- CurData[order(CurData$x),]
          plot(x = NA, y = NA, xlim = c(xLwr, xUpr), ylim = c(yLwr, yUpr),
               main = "", xlab = "", ylab = "", axes = FALSE)
          lines(x = CurData$x, y = CurData$Mismatch, lwd = LineWidth)
          lines(x = CurData$x, y = CurData$lwr, lty = 3)
          lines(x = CurData$x, y = CurData$upr, lty = 3)
          # Add in the Zopt line
          lines(x = CurData$x, y = CurData$Zopt, col = "grey")
          box()
     }
     dev.off()
}
# Create the gif in the main directory and then remove the temporary directory
GifName <- "ResultFigures/Animations/Equilibrium/MuNiche.gif"
SysCommand1 <- paste("convert -delay 50 ", TempDirPath, "/gen*.png ", GifName, sep = "")
SysCommand2 <- paste("rm -r ", TempDirPath, sep = "")
system(SysCommand1)
system(SysCommand2)
rm(MismatchData)

# --------------------------------------------- 3) genetic variance in dispersal
load("SimData/EquilibriumDispGenVarData.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- Inf
xUpr <- -Inf
yLwr <- Inf
yUpr <- -Inf
for(p in 1:9){
     CurData <- DispGenVarData[[p]]
     xLwr <- min(xLwr, min(CurData$x))
     xUpr <- max(xUpr, max(CurData$x))
     yLwr <- min(yLwr, min(CurData$lwr))
     yUpr <- max(yUpr, max(CurData$upr))
}

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
     par(mfrow = c(3,3), oma = OuterMar, mar = InnerMar)
     for(i in 1:9){
          CurData <- subset(DispGenVarData[[i]], g == gen)
          CurData <- CurData[order(CurData$x),]
          plot(x = NA, y = NA, xlim = c(xLwr, xUpr), ylim = c(yLwr, yUpr),
               main = "", xlab = "", ylab = "", axes = FALSE)
          lines(x = CurData$x, y = CurData$GenVar, lwd = LineWidth)
          lines(x = CurData$x, y = CurData$lwr, lty = 3)
          lines(x = CurData$x, y = CurData$upr, lty = 3)
          box()
     }
     dev.off()
}
# Create the gif in the main directory and then remove the temporary directory
GifName <- "ResultFigures/Animations/Equilibrium/DispGenVar.gif"
SysCommand1 <- paste("convert -delay 50 ", TempDirPath, "/gen*.png ", GifName, sep = "")
SysCommand2 <- paste("rm -r ", TempDirPath, sep = "")
system(SysCommand1)
system(SysCommand2)
rm(DispGenVarData)

# ------------------------------------------- 3) genetic variance in niche trait
load("SimData/EquilibriumFitGenVarData.rdata")
# Figure out appropriate values for xlim and ylim
xLwr <- Inf
xUpr <- -Inf
yLwr <- Inf
yUpr <- -Inf
for(p in 1:9){
     CurData <- FitGenVarData[[p]]
     xLwr <- min(xLwr, min(CurData$x))
     xUpr <- max(xUpr, max(CurData$x))
     yLwr <- min(yLwr, min(CurData$lwr))
     yUpr <- max(yUpr, max(CurData$upr))
}

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
     par(mfrow = c(3,3), oma = OuterMar, mar = InnerMar)
     for(i in 1:9){
          CurData <- subset(FitGenVarData[[i]], g == gen)
          CurData <- CurData[order(CurData$x),]
          plot(x = NA, y = NA, xlim = c(xLwr, xUpr), ylim = c(yLwr, yUpr),
               main = "", xlab = "", ylab = "", axes = FALSE)
          lines(x = CurData$x, y = CurData$GenVar, lwd = LineWidth)
          lines(x = CurData$x, y = CurData$lwr, lty = 3)
          lines(x = CurData$x, y = CurData$upr, lty = 3)
          box()
     }
     dev.off()
}
# Create the gif in the main directory and then remove the temporary directory
GifName <- "ResultFigures/Animations/Equilibrium/FitGenVar.gif"
SysCommand1 <- paste("convert -delay 50 ", TempDirPath, "/gen*.png ", GifName, sep = "")
SysCommand2 <- paste("rm -r ", TempDirPath, sep = "")
system(SysCommand1)
system(SysCommand2)
rm(FitGenVarData)
