# This script will finally create the graphs from the abundance data.

# Set the working directory and declare which paramater combination to work with
setwd("~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/")
RangeParams <- 1

# Load in the abundance data
InFile <- paste("SimData/Params", RangeParams, "ExtractedAbunds.rdata", sep = "")
load(InFile)

# First make a graph of the SectorMeans and SectorSigma through time
library(plot3D)

pdf(file = "SectorAbund.pdf", width = 7, height = 3, onefile = FALSE, paper = "special")
     par(mfrow = c(1, 2))
     # First the mean abundances
     image2D(z = SectorMeans, xaxt = "n", yaxt = "n", xlab = "Sector", ylab = "Time",
             main = "Mean Abundance")
     axis(1)
     axis(2)
     
     # Now the standard deviations
     image2D(z = SectorSigma, xaxt = "n", yaxt = "n", xlab = "Sector", ylab = "Time",
             main = "Standard Deviation of Abundances")
     axis(1)
     axis(2)
dev.off()

# Now make the animated gifs of abundance in each patch through time
nTime <- dim(MeanAbunds)[3]
TimeSeq <- seq(1, nTime, by = 10)

# First make a convenience function and establish some constant color pallettes to use
MeanRange <- range(MeanAbunds)
SigmaRange <- range(SigmaAbunds)
MeanCols <- data.frame(Index = seq(MeanRange[1], MeanRange[2], length.out = 10000),
                       Cols = jet.col(10000))
SigmaCols <- data.frame(Index = seq(SigmaRange[1], SigmaRange[2], length.out = 10000),
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
#    combined into a single gif
for(t in TimeSeq){
     PlotFile <- paste("TempGifAssembly/AbundGen", t, ".pdf", sep = "")
     CurMeanAbunds <- MeanAbunds[,,t]
     CurSigmaAbunds <- SigmaAbunds[,,t]
     pdf(file = PlotFile, width = 7, height = 3, onefile = FALSE, paper = "special")
          par(mfrow = c(1, 2))
          # First the mean abundances
          ColRange <- FindRange(minimum = min(CurMeanAbunds), maximum = max(CurMeanAbunds),
                                sequence = MeanCols$Index)
          image2D(z = CurMeanAbunds, xaxt = "n", yaxt = "n", xlab = "Width", ylab = "Length",
                  main = "Mean Abundance", col = MeanCols$Cols[ColRange])
          axis(1)
          axis(2)
     
          # Now the standard deviations
          ColRange <- FindRange(minimum = min(CurSigmaAbunds), maximum = max(CurSigmaAbunds),
                                sequence = SigmaCols$Index)
          image2D(z = CurSigmaAbunds, xaxt = "n", yaxt = "n", xlab = "Width", ylab = "Length",
                  main = "Standard Deviation of Abundances", col = SigmaCols$Cols[ColRange])
          axis(1)
          axis(2)
          GenText <- paste("Generation ", t, sep = "")
          text(x = 10, y = 100, labels = GenText, xpd = NA)
     dev.off()
}

# Finally, make a complete gif and then clean up all the single images
GifName <- paste("AbundGif", RangeParams, ".gif", sep = "")
SysCommand <- paste("convert -delay 25 TempGifAssembly/AbundGen* ", GifName, sep = "")
system(SysCommand)


# If the gif is satisfactory, then run this line of code to clean up the TempGifAssembly
#    folder before trying this for another set of range parameters
system("rm TempGifAssembly/AbundGen*")
