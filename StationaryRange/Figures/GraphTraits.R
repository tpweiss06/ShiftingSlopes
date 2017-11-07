# This script will finally create the graphs from the abundance data.

# Set the working directory and declare which paramater combination to work with
setwd("~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/")
RangeParams <- 1

# Load in the abundance data
InFile <- paste("SimData/Params", RangeParams, "/ExtractedTraits.rdata", sep = "")
load(InFile)

# First make a graph of the SectorMeans and SectorSigma through time
library(plot3D)
LenSeq <- 450:550
LocLabels <- seq(-50, 50, by = 20)
TimeSeq <- 1:1000 
TimeLabels <- seq(0, 1000, by = 200)
TextSize <- 1.5
AxisSize <- 1.1
FitFile <- paste("Figures/Params", RangeParams, "/SectorFit.pdf", sep = "")
DispFile <- paste("Figures/Params", RangeParams, "/SectorDisp.pdf", sep = "")
pdf(file = FitFile, width = 8, height = 5, onefile = FALSE, paper = "special")
     par(mfrow = c(1, 2), mar = c(2.5, 2.75, 2, 2.75) + 0.05, oma = c(1,1,1,0))
     # First the mean fitness
     image2D(z = SectorMuFit[LenSeq,TimeSeq], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             main = "")
     mtext("Mean", side = 3, cex = TextSize)
     axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
     axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
     axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     
     # Now the standard deviations
     image2D(z = SectorSigFit[LenSeq,TimeSeq], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             main = "")
     mtext("Standard deviation", side = 3, cex = TextSize)
     axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
     axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
     axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     
     # Axis labels
     mtext("Sector", side = 1, outer = TRUE, line = -0.5, cex = TextSize)
     mtext("Generation", side = 2, outer = TRUE, line = -0.25, cex = TextSize)
     mtext("Fitness", side = 3, outer = TRUE, line = -0.5, cex = TextSize)
dev.off()

pdf(file = DispFile, width = 8, height = 5, onefile = FALSE, paper = "special")
     par(mfrow = c(1, 2), mar = c(2.5, 2.75, 2, 2.75) + 0.05, oma = c(1,1,1,0))
     # First the mean fitness
     image2D(z = SectorMuDisp[LenSeq,TimeSeq], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             main = "")
     mtext("Mean", side = 3, cex = TextSize)
     axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
     axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
     axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)

     # Now the standard deviations
     image2D(z = SectorSigDisp[LenSeq,TimeSeq], xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             main = "")
     mtext("Standard deviation", side = 3, cex = TextSize)
     axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
     axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
     axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)

     # Axis labels
     mtext("Sector", side = 1, outer = TRUE, line = -0.5, cex = TextSize)
     mtext("Generation", side = 2, outer = TRUE, line = -0.25, cex = TextSize)
     mtext("Dispersal", side = 3, outer = TRUE, line = -0.5, cex = TextSize)
dev.off()

# Now make the animated gifs of abundance in each patch through time
# First make a convenience function and establish some constant color pallettes to use
MuFitRange <- range(MeanMuFit, na.rm = TRUE)
SigFitRange <- range(MeanSigFit, na.rm = TRUE)
MuDispRange <- range(MeanMuDisp, na.rm = TRUE)
SigDispRange <- range(MeanSigDisp, na.rm = TRUE)
HeatMapCols <- colorRampPalette(c("blue", "yellow", "red"))
MuFitCols <- data.frame(Index = seq(MuFitRange[1], MuFitRange[2], length.out = 10000),
                        Cols = HeatMapCols(10000))
SigFitCols <- data.frame(Index = seq(SigFitRange[1], SigFitRange[2], length.out = 10000),
                         Cols = HeatMapCols(10000))
MuDispCols <- data.frame(Index = seq(MuDispRange[1], MuDispRange[2], length.out = 10000),
                         Cols = HeatMapCols(10000))
SigDispCols <- data.frame(Index = seq(SigDispRange[1], SigDispRange[2], length.out = 10000),
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

# Create a new TimeSeq vector with an entry for every 20 generations
TimeSeq <- seq(1, 1000, by = 20)
WidthLabels <- seq(1, 10, by = 1)

# Now use those objects to create the single generation images that will be
#    combined into a single gif for fitness
for(t in TimeSeq){
     if(t < 10){
          PlotFile <- paste("Figures/TempGifAssembly/FitGen00", t, ".png", sep = "")
     } else if(t < 100){
          PlotFile <- paste("Figures/TempGifAssembly/FitGen0", t, ".png", sep = "")
     } else{
          PlotFile <- paste("Figures/TempGifAssembly/FitGen", t, ".png", sep = "")
     }
     CurMuFit <- MeanMuFit[,LenSeq,t]
     CurSigFit <- MeanSigFit[,LenSeq,t]
     png(filename = PlotFile, width = 8, height = 5, units = "in", pointsize = 12, res = 500)
          par(mfrow = c(1, 2), mar = c(2.5, 3.25, 2, 3.25) + 0.05, oma = c(1,1.25,1,0.5))
          # First the means
          ColRange <- FindRange(minimum = min(CurMuFit, na.rm = TRUE), maximum = max(CurMuFit, na.rm = TRUE),
                                sequence = MuFitCols$Index)
          image2D(z = CurMuFit, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = "", col = as.character(MuFitCols$Cols[ColRange]), colkey = FALSE)
          # Add the color key
          colkey(col = as.character(MuFitCols$Cols), clim = MuFitRange, 
                 cex.axis = AxisSize, add = TRUE)
          mtext("Mean", side = 3, cex = TextSize)
          axis(1, at = seq(0.05, 0.95, by = 0.1), labels = WidthLabels, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.2), labels = LocLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          
          # Now the standard deviations
          ColRange <- FindRange(minimum = min(CurSigFit, na.rm = TRUE), maximum = max(CurSigFit, na.rm = TRUE),
                                sequence = SigFitCols$Index)
          image2D(z = CurSigFit, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = "", col = as.character(SigFitCols$Cols[ColRange]), colkey = FALSE)
          # Add the color key
          colkey(col = as.character(SigFitCols$Cols), clim = SigFitRange, 
                 cex.axis = AxisSize, add = TRUE)
          mtext("Standard deviation", side = 3, cex = TextSize)
          axis(1, at = seq(0.05, 0.95, by = 0.1), labels = WidthLabels, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          
          # Axis labels
          mtext("Width", side = 1, outer = TRUE, line = -0.5, cex = TextSize)
          mtext("Length", side = 2, outer = TRUE, line = -0.25, cex = TextSize)
          mtext("Fitness", side = 3, outer = TRUE, line = -0.5, cex = TextSize)
          
          # Generation text
          GenText <- paste("Generation ", t, sep = "")
          text(x = 1.1, y = 1.1, labels = GenText, xpd = NA)
     dev.off()
}

# Finally, make a complete gif and then clean up all the single images
GifName <- paste("~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/Figures/Params", RangeParams, "/FitGif.gif", sep = "")
SysCommand <- paste("convert -delay 50 ~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/Figures/TempGifAssembly/FitGen* ", GifName, sep = "")
system(SysCommand)

# If the gif is satisfactory, then run this line of code to clean up the TempGifAssembly
#    folder before trying this for another set of range parameters
system("rm ~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/Figures/TempGifAssembly/FitGen*")


# Now use those objects to create the single generation images that will be
#    combined into a single gif for dispersal
for(t in TimeSeq){
     if(t < 10){
          PlotFile <- paste("Figures/TempGifAssembly/DispGen00", t, ".png", sep = "")
     } else if(t < 100){
          PlotFile <- paste("Figures/TempGifAssembly/DispGen0", t, ".png", sep = "")
     } else{
          PlotFile <- paste("Figures/TempGifAssembly/DispGen", t, ".png", sep = "")
     }
     CurMuDisp <- MeanMuDisp[,LenSeq,t]
     CurSigDisp <- MeanSigDisp[,LenSeq,t]
     png(filename = PlotFile, width = 8, height = 5, units = "in", pointsize = 12, res = 500)
          par(mfrow = c(1, 2), mar = c(2.5, 3.25, 2, 3.25) + 0.05, oma = c(1,1.25,1,0.5))
          # First the means
          ColRange <- FindRange(minimum = min(CurMuDisp, na.rm = TRUE), maximum = max(CurMuDisp, na.rm = TRUE),
                                sequence = MuDispCols$Index)
          image2D(z = CurMuDisp, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = "", col = as.character(MuDispCols$Cols[ColRange]), colkey = FALSE)
          # Add the color key
          colkey(col = as.character(MuDispCols$Cols), clim = MuDispRange, 
                 cex.axis = AxisSize, add = TRUE)
          mtext("Mean", side = 3, cex = TextSize)
          axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
          axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     
          # Now the standard deviations
          ColRange <- FindRange(minimum = min(CurSigDisp, na.rm = TRUE), maximum = max(CurSigDisp, na.rm = TRUE),
                                sequence = SigDispCols$Index)
          image2D(z = CurSigDisp, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = "", col = as.character(SigDispCols$Cols[ColRange]), colkey = FALSE)
          # Add the color key
          colkey(col = as.character(SigDispCols$Cols), clim = SigDispRange, 
                 cex.axis = AxisSize, add = TRUE)
          mtext("Standard deviation", side = 3, cex = TextSize)
          axis(1, at = seq(0, 1, by = 0.2), labels = LocLabels, cex.axis = AxisSize)
          axis(1, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          axis(2, at = seq(0, 1, by = 0.2), labels = TimeLabels, las = 1, cex.axis = AxisSize)
          axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     
          # Axis labels
          mtext("Width", side = 1, outer = TRUE, line = -0.5, cex = TextSize)
          mtext("Length", side = 2, outer = TRUE, line = -0.5, cex = TextSize)
          mtext("Dispersal", side = 3, outer = TRUE, line = -0.5, cex = TextSize)
     
          # Generation text
          GenText <- paste("Generation ", t, sep = "")
          text(x = 1.1, y = 1.1, labels = GenText, xpd = NA)
     dev.off()
}

# Finally, make a complete gif and then clean up all the single images
GifName <- paste("~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/Figures/Params", RangeParams, "/DispGif.gif", sep = "")
SysCommand <- paste("convert -delay 50 ~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/Figures/TempGifAssembly/DispGen* ", GifName, sep = "")
system(SysCommand)

# If the gif is satisfactory, then run this line of code to clean up the TempGifAssembly
#    folder before trying this for another set of range parameters
system("rm ~/Desktop/RangeShifts/ShiftingSlopes/StationaryRange/Figures/TempGifAssembly/DispGen*")



