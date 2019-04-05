# Plot the var data. Doesn't have to be pretty, just has to let me see a visualization
#    and then get some key values

load("~/Desktop/RangeShifts/ShiftingSlopesOther/SimData/VarData.rdata")
# Determine the x range to use for all graphs
xMin <- 0
xMax <- 0
for(p in 1:9){
     xRange <- range(VarList[[p]]$x, na.rm = TRUE)
     xMin <- min(c(xRange[1], xMin))
     xMax <- max(c(xRange[2], xMax))
}
xRange <- c(xMin, xMax)

# First make the graph for variation in abundances
quartz(width = 12, height = 9)
par(mfrow = c(3,3))
for(p in 1:9){
     plot(x = VarList[[p]]$x, y = VarList[[p]]$abund, xlim = xRange, xlab = "Space", ylab = "Var(abund)")
}

# Now make the graph for variation in mean niche trait
quartz(width = 12, height = 9)
par(mfrow = c(3,3))
for(p in 1:9){
     plot(x = VarList[[p]]$x, y = VarList[[p]]$muFit, xlim = xRange, xlab = "Space", ylab = "Var(niche)")
}

# Now make the graph for variation in mean dispersal
quartz(width = 12, height = 9)
par(mfrow = c(3,3))
for(p in 1:9){
     plot(x = VarList[[p]]$x, y = VarList[[p]]$muDisp, xlim = xRange, xlab = "Space", ylab = "Var(disp)")
}

# Now make the graph for variation in variance in fitness of patches
quartz(width = 12, height = 9)
par(mfrow = c(3,3))
for(p in 1:9){
     plot(x = VarList[[p]]$x, y = VarList[[p]]$sigmaFitPhen, xlim = xRange, xlab = "Space", ylab = "Patch Var(fit)")
}


# Look into the maximum standard deviations across the y dimension
MaxAbundSD <- rep(NA, 9)
MaxNicheSD <- rep(NA, 9)
MaxDispSD <- rep(NA, 9)
for(p in 1:9){
     MaxAbundSD[p] <- (max(VarList[[p]]$abund, na.rm = TRUE))
     MaxNicheSD[p] <- (max(VarList[[p]]$muFit, na.rm = TRUE))
     MaxDispSD[p] <- (max(VarList[[p]]$muDisp, na.rm = TRUE))
}
MaxAbundSD
MaxNicheSD
MaxDispSD


#  Plot the average dispersal CV across simulations for all x points
quartz(width = 12, height = 9)
par(mfrow = c(3,3))
for(p in 1:9){
     xSeq <- unique(VarList[[p]]$x)
     DispCVs <- rep(NA, length(xSeq))
     for(i in 1:length(xSeq)){
          CurData <- subset(VarList[[p]], x == xSeq[i])
          DispCVs[i] <- mean(CurData$muDisp, na.rm = TRUE)
     }
     plot(x = xSeq, y = DispCVs, xlim = xRange, xlab = "Space", ylab = "CV disp",
          ylim = c(0, 0.03))
     abline(h = 0.01)
     abline(h = 0.02)
     abline(h = 0.03)
}


#  Plot the average dispersal CV across simulations for all x points
quartz(width = 12, height = 9)
par(mfrow = c(3,3))
for(p in 1:9){
     xSeq <- unique(VarList[[p]]$x)
     FitCVs <- rep(NA, length(xSeq))
     for(i in 1:length(xSeq)){
          CurData <- subset(VarList[[p]], x == xSeq[i])
          FitCVs[i] <- mean(CurData$muFit, na.rm = TRUE)
     }
     plot(x = xSeq, y = FitCVs, xlim = xRange, xlab = "Space", ylab = "CV fit",
          ylim = c(-0.5, 0.5))
     abline(h = 0)
     abline(h = 0.25, lty = 2)
     abline(h = -0.25, lty = 2)
}

