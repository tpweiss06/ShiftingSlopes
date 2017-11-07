# This script will finally create the graphs from the abundance data.

# Set the working directory and declare which paramater combination to work with
setwd("~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/")
RangeParams <- 2

# Load in the abundance data
InFile <- paste("SimData/Params", RangeParams, "/ExtractedAbunds.rdata", sep = "")
load(InFile)

# First make a graph of the SectorMeans and SectorSigma through time
library(plot3D)
LenSeq <- 450:550
LocLabels <- seq(-50, 50, by = 20)
TimeSeq <- 1:50
TimeLabels <- seq(0, 50, by = 10)
TextSize <- 1.5
AxisSize <- 1.1
FigFile <- paste("Figures/Params", RangeParams, "/SectorAbund.pdf", sep = "")
pdf(file = FigFile, width = 8, height = 5, onefile = FALSE, paper = "special")
     par(mfrow = c(1, 2), mar = c(2.5, 2.75, 2, 2.75) + 0.05, oma = c(1,1,1,0))
     # First the mean abundances
     image2D(z = SectorMeans[LenSeq,TimeSeq], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             main = "")
     mtext("Mean", side = 3, cex = TextSize)
     axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
     axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
     axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     
     # Now the standard deviations
     image2D(z = SectorSigma[LenSeq,TimeSeq], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             main = "")
     mtext("Standard deviation", side = 3, cex = TextSize)
     axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
     axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
     axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     
     # Axis labels
     mtext("Sector", side = 1, outer = TRUE, line = -0.5, cex = TextSize)
     mtext("Generation", side = 2, outer = TRUE, line = -0.5, cex = TextSize)
     mtext("Abundance", side = 3, outer = TRUE, line = -0.5, cex = TextSize)
dev.off()

# Now make the animated gifs of abundance in each patch through time
# First make a convenience function and establish some constant color pallettes to use
MeanRange <- range(MeanAbunds)
SigmaRange <- range(SigmaAbunds)
HeatMapCols <- colorRampPalette(c("blue", "yellow", "red"))
MeanCols <- data.frame(Index = seq(MeanRange[1], MeanRange[2], length.out = 10000),
                       Cols = HeatMapCols(10000))
SigmaCols <- data.frame(Index = seq(SigmaRange[1], SigmaRange[2], length.out = 10000),
                       Cols = HeatMapCols(10000))

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
WidthLabels <- seq(1, 10, by = 1)

for(t in TimeSeq){
     if(t < 10){
          PlotFile <- paste("Figures/TempGifAssembly/AbundGen0", t, ".png", sep = "")
     } else{
          PlotFile <- paste("Figures/TempGifAssembly/AbundGen", t, ".png", sep = "")
     }
     CurMeanAbunds <- MeanAbunds[,LenSeq,t]
     CurSigmaAbunds <- SigmaAbunds[,LenSeq,t]
     png(filename = PlotFile, width = 8, height = 5, units = "in", pointsize = 12, res = 500)
          par(mfrow = c(1, 2), mar = c(2.5, 3, 2, 3) + 0.05, oma = c(1,1,1,0.5))
          # First the mean abundances
          ColRange <- FindRange(minimum = min(CurMeanAbunds), maximum = max(CurMeanAbunds),
                                sequence = MeanCols$Index)
          image2D(z = CurMeanAbunds, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = "", col = as.character(MeanCols$Cols[ColRange]), colkey = FALSE)
          # Add the color key
          colkey(col = as.character(MeanCols$Cols), clim = MeanRange, 
                 cex.axis = AxisSize, add = TRUE)
          mtext("Mean", side = 3, cex = TextSize)
          axis(1, at = seq(0.05, 0.95, by = 0.1), labels = WidthLabels, cex.axis = AxisSize)          
          axis(2, at = seq(0, 1, by = 0.2), labels = LocLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     
          # Now the standard deviations
          ColRange <- FindRange(minimum = min(CurSigmaAbunds), maximum = max(CurSigmaAbunds),
                                sequence = SigmaCols$Index)
          image2D(z = CurSigmaAbunds, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = "", col = as.character(SigmaCols$Cols[ColRange]),
                  colkey = FALSE)
          # Add the color key
          colkey(col = as.character(SigmaCols$Cols), clim = SigmaRange, 
                 cex.axis = AxisSize, add = TRUE)
          mtext("Standard deviation", side = 3, cex = TextSize)
          axis(1, at = seq(0.05, 0.95, by = 0.1), labels = WidthLabels, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.2), labels = LocLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          
          # Axis labels
          mtext("Width", side = 1, outer = TRUE, line = -0.5, cex = TextSize)
          mtext("Length", side = 2, outer = TRUE, line = -0.5, cex = TextSize)
          mtext("Abundance", side = 3, outer = TRUE, line = -0.5, cex = TextSize)
          
          # Generation text
          GenText <- paste("Generation ", t, sep = "")
          text(x = 1.1, y = 1.1, labels = GenText, xpd = NA)
     dev.off()
}

# Finally, make a complete gif and then clean up all the single images
GifName <- paste("~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/Figures/Params", RangeParams, "/AbundGif.gif", sep = "")
SysCommand <- paste("convert -delay 50 ~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/Figures/TempGifAssembly/AbundGen* ", GifName, sep = "")
system(SysCommand)

# If the gif is satisfactory, then run this line of code to clean up the TempGifAssembly
#    folder before trying this for another set of range parameters
system("rm ~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/Figures/TempGifAssembly/AbundGen*")
