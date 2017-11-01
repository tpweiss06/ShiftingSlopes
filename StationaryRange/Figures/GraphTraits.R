# This script will finally create the graphs from the abundance data.

# Set the working directory and declare which paramater combination to work with
setwd("~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/")
RangeParams <- 1

# Load in the abundance data
InFile <- paste("SimData/Params", RangeParams, "ExtractedTraits.rdata", sep = "")
load(InFile)

# First make a graph of the SectorMeans and SectorSigma through time
library(plot3D)

pdf(file = "SectorFit.pdf", width = 7, height = 3, onefile = FALSE, paper = "special")
     par(mfrow = c(1, 2))
     # First the mean abundances
     image2D(z = SectorMuFit, xaxt = "n", yaxt = "n", xlab = "Sector", ylab = "Time",
             main = "Mean Fitness")
     axis(1)
     axis(2)

     # Now the standard deviations
     image2D(z = SectorSigFit, xaxt = "n", yaxt = "n", xlab = "Sector", ylab = "Time",
             main = "Standard Deviation of Fitness")
     axis(1)
     axis(2)
dev.off()

pdf(file = "SectorDisp.pdf", width = 7, height = 3, onefile = FALSE, paper = "special")
     par(mfrow = c(1, 2))
     # First the means
     image2D(z = SectorMuDisp, xaxt = "n", yaxt = "n", xlab = "Sector", ylab = "Time",
             main = "Mean Dispersal")
     axis(1)
     axis(2)

     # Now the standard deviations
     image2D(z = SectorSigDisp, xaxt = "n", yaxt = "n", xlab = "Sector", ylab = "Time",
             main = "Standard Deviation of Dispersal")
     axis(1)
     axis(2)
dev.off()

# Now make the animated gifs of abundance in each patch through time
nTime <- dim(MeanMuFit)[3]
TimeSeq <- seq(1, nTime, by = 10)

# First make a convenience function and establish some constant color pallettes to use
MuFitRange <- range(MeanMuFit)
SigFitRange <- range(MeanSigFit)
MuDispRange <- range(MeanMuDisp)
SigDispRange <- range(MeanSigDisp)

MuFitCols <- data.frame(Index = seq(MuFitRange[1], MuFitRange[2], length.out = 10000),
                        Cols = jet.col(10000))
SigFitCols <- data.frame(Index = seq(SigFitRange[1], SigFitRange[2], length.out = 10000),
                         Cols = jet.col(10000))
MuDispCols <- data.frame(Index = seq(MuDispRange[1], MuDispRange[2], length.out = 10000),
                         Cols = jet.col(10000))
SigDispCols <- data.frame(Index = seq(SigDispRange[1], SigDispRange[2], length.out = 10000),
                          Cols = jet.col(10000))

# Create a function to find the indices corresponding to the subset of the a sequence
#	which is defined as the shortest possible sequence including both the given
#	minimum and maximum
FindRange <- function(minimum, maximum, sequence){
     LessThan <- ifelse(sequence <= minimum, 1, 0)
     MoreThan <- ifelse(sequence >= maximum, 1, 0)
     MinIndex <- sum(LessThan)
     MaxIndex <- length(sequence) - sum(MoreThan)
     return(MinIndex:MaxIndex)
}

# Now use those objects to create the single generation images that will be
#    combined into a single gif for fitness
for(t in TimeSeq){
     PlotFile <- paste("TempGifAssembly/FitGen", t, ".pdf", sep = "")
     CurMuFit <- MeanMuFit[,,t]
     CurSigFit <- MeanSigFit[,,t]
     pdf(file = PlotFile, width = 7, height = 3, onefile = FALSE, paper = "special")
     par(mfrow = c(1, 2))
     # First the means
     ColRange <- FindRange(minimum = min(CurMuFit), maximum = max(CurMuFit),
                           sequence = MuFitCols$Index)
     image2D(z = CurMuFit, xaxt = "n", yaxt = "n", xlab = "Width", ylab = "Length",
             main = "Mean Fitness", col = MuFitCols$Cols[ColRange])
     axis(1)
     axis(2)
     
     # Now the standard deviations
     ColRange <- FindRange(minimum = min(CurSigFit), maximum = max(CurSigFit),
                           sequence = SigFitCols$Index)
     image2D(z = CurSigFit, xaxt = "n", yaxt = "n", xlab = "Width", ylab = "Length",
             main = "Standard Deviation of Fitness", col = SigFitCols$Cols[ColRange])
     axis(1)
     axis(2)
     GenText <- paste("Generation ", t, sep = "")
     text(x = 10, y = 100, labels = GenText, xpd = NA)
     dev.off()
}

# Finally, make a complete gif and then clean up all the single images
GifName <- paste("FitGif", RangeParams, ".gif", sep = "")
SysCommand <- paste("convert -delay 25 TempGifAssembly/FitGen* ", GifName, sep = "")
system(SysCommand)


# If the gif is satisfactory, then run this line of code to clean up the TempGifAssembly
#    folder before trying this for another set of range parameters
system("rm TempGifAssembly/FitGen*")


# Now use those objects to create the single generation images that will be
#    combined into a single gif for dispersal
for(t in TimeSeq){
     PlotFile <- paste("TempGifAssembly/DispGen", t, ".pdf", sep = "")
     CurMuFit <- MeanMuDisp[,,t]
     CurSigFit <- MeanSigDisp[,,t]
     pdf(file = PlotFile, width = 7, height = 3, onefile = FALSE, paper = "special")
          par(mfrow = c(1, 2))
          # First the means
          ColRange <- FindRange(minimum = min(CurMuDisp), maximum = max(CurMuDisp),
                                sequence = MuDispCols$Index)
          image2D(z = CurMuDisp, xaxt = "n", yaxt = "n", xlab = "Width", ylab = "Length",
                  main = "Mean Dispersal", col = MuDispCols$Cols[ColRange])
          axis(1)
          axis(2)
     
          # Now the standard deviations
          ColRange <- FindRange(minimum = min(CurSigDisp), maximum = max(CurSigDisp),
                                sequence = SigDispCols$Index)
          image2D(z = CurSigDisp, xaxt = "n", yaxt = "n", xlab = "Width", ylab = "Length",
                  main = "Standard Deviation of Dispersal", col = SigDispCols$Cols[ColRange])
          axis(1)
          axis(2)
          GenText <- paste("Generation ", t, sep = "")
          text(x = 10, y = 100, labels = GenText, xpd = NA)
     dev.off()
}

# Finally, make a complete gif and then clean up all the single images
GifName <- paste("DispGif", RangeParams, ".gif", sep = "")
SysCommand <- paste("convert -delay 25 TempGifAssembly/DispGen* ", GifName, sep = "")
system(SysCommand)


# If the gif is satisfactory, then run this line of code to clean up the TempGifAssembly
#    folder before trying this for another set of range parameters
system("rm TempGifAssembly/DispGen*")



