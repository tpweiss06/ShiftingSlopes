# Generate a figure demonstrating the idea of spatial population structure and
#    local adaptation.

# First create some of the colors I will need for the graph
LocalAdaptCols <- colorRampPalette(c("green", "blue", "red"))
EnvGradCols <- colorRampPalette(c("black", "white", "black"))

# Now create data frames to hold the x and y coordinates for graphing
NoStruct <- data.frame(x = rep(NA, 78), y = rep(NA, 78))
SpatStruct <- NoStruct

NoStruct$x <- rep(seq(-6, 6, by = 1), each = 6)
NoStruct$y <- rep(seq(1, 6, by = 1), 13)
SpatStruct$x <- c(-8.5, rep(-7.5, 2), rep(-6.5, 3), rep(-5.5, 4), rep(-4.5, 5),
                  rep(seq(-3.5, 3.5, by = 1), each = 6), 
                  rep(4.5, 5), rep(5.5, 4), rep(6.5, 3), rep(7.5, 2), 8.5)
SpatStruct$y <- c(1, 1:2, 1:3, 1:4, 1:5, rep(1:6, 8), 1:5, 1:4, 1:3, 1:2, 1)


AdaptCols <- LocalAdaptCols(18)
AdaptColSeq <- NULL
for(i in 1:5){
     AdaptColSeq <- c(AdaptColSeq, rep(AdaptCols[i], i))
}
for(i in 6:13){
     AdaptColSeq <- c(AdaptColSeq, rep(AdaptCols[i], 6))
}
for(i in 14:18){
     AdaptColSeq <- c(AdaptColSeq, rep(AdaptCols[i], i - (9 + 2*(14-i))))
}

pch_seq <- c(rep(15, 200), rep(0, 600), rep(15, 200))

pdf(file = "~/Desktop/RangeShifts/ShiftingSlopes/SpatialPopulationStructure.pdf", 
    width = 4, height = 3, onefile = FALSE, paper = "special")
     par(mfrow = c(2,1), mar = c(0,3,0,3))
     plot(SpatStruct$x, SpatStruct$y, xlim = c(-10, 10), ylim = c(-2, 8), cex = 1,
          col = AdaptColSeq, pch = 19, axes = FALSE, xlab = "", ylab = "")
     polygon(x = c(-15, -15, 15, 15), y = c(0, -1.5, -1.5, 0), xpd = NA)
     nSlices <- 1000
     PolCols <- EnvGradCols(nSlices)
     xSlices <- seq(-15, 15, length.out = nSlices + 1)
     for(i in 1:nSlices){
          polygon(x = c(xSlices[i], xSlices[i], xSlices[i+1], xSlices[i+1]),
                  y = c(0, -1.5, -1.5, 0), col = PolCols[i], border = NA, xpd = NA)
     }
     text("(a)", x = -14, y = 7.5, cex = 0.75, xpd = NA)
     
     plot(NoStruct$x, NoStruct$y, xlim = c(-10, 10), ylim = c(-2, 8), cex = 1,
          col = sample(AdaptColSeq, size = length(AdaptColSeq), replace = FALSE), 
          pch = 19, axes = FALSE, xlab = "", ylab = "")
     polygon(x = c(-15, -15, -6.5, -6.5), y = c(0, -1.5, -1.5, 0), xpd = NA, col = "black")
     polygon(x = c(-6.5, -6.5, 6.5, 6.5), y = c(0, -1.5, -1.5, 0), xpd = NA, col = "white")
     polygon(x = c(6.5, 6.5, 15, 15), y = c(0, -1.5, -1.5, 0), xpd = NA, col = "black")
     text("(b)", x = -14, y = 7.5, cex = 0.75, xpd = NA)
dev.off()

pdf(file = "~/Desktop/RangeShifts/ShiftingSlopesOther/SchematicFigures/LocalAdaptGrad.pdf", 
    width = 4, height = 3, onefile = FALSE, paper = "special")
     par(mfrow = c(2,1), mar = c(0,3,0,3))
     plot(SpatStruct$x, SpatStruct$y, xlim = c(-10, 10), ylim = c(-2, 8), cex = 1,
          col = AdaptColSeq, pch = 19, axes = FALSE, xlab = "", ylab = "")
     
     plot(SpatStruct$x, SpatStruct$y, xlim = c(-10, 10), ylim = c(-2, 8), cex = 1,
          col = sample(AdaptColSeq, size = length(AdaptColSeq), replace = FALSE), 
          pch = 19, axes = FALSE, xlab = "", ylab = "")
dev.off()

pdf(file = "~/Desktop/RangeShifts/ShiftingSlopesOther/SchematicFigures/EdgeGrad.pdf", 
    width = 4, height = 3, onefile = FALSE, paper = "special")
     par(mfrow = c(2,1), mar = c(0,3,0,3))
     plot(SpatStruct$x, SpatStruct$y, xlim = c(-10, 10), ylim = c(-2, 8), cex = 1,
          col = "black", pch = 19, axes = FALSE, xlab = "", ylab = "")

     plot(NoStruct$x, NoStruct$y, xlim = c(-10, 10), ylim = c(-2, 8), cex = 1,
          col = "black", pch = 19, axes = FALSE, xlab = "", ylab = "")
dev.off()