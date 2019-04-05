# This script will create the simplified composite graph incorporating dispersal
#    evolution and initial dispersal

setwd("~/Desktop/RangeShifts/ShiftingSlopesOther/ResultFigures/")
load("SimData/InitDispData.rdata")
load("SimData/DispEvolData.rdata")
CurSpeed <- 2
Params <- c(1, 9)

AllSims <- vector(length = 2, mode = "list")
ExtantSims <- vector(length = 2, mode = "list")
NumExtant <- rep(0, 2)
NumSims <- 200
xMin <- 0
xMax <- 0
for(p in 1:2){
     CurParam <- Params[p]
     AllVals <- NULL
     ExtantVals <- NULL
     for(i in 1:NumSims){
          AllVals <- c(AllVals, DispInit[[CurParam]][[i]]$ExpDists)
          
          if(DispInit[[CurParam]][[i]]$Ext[CurSpeed] == 1){
               NumExtant[p] <- NumExtant[p] + 1
               ExtantVals <- c(ExtantVals, DispInit[[CurParam]][[i]]$ExpDists)
          }
     }
     AllSims[[CurParam]] <- log(AllVals, base = 10)
     if(!is.null(ExtantVals)){
          ExtantSims[[CurParam]] <- log(ExtantVals, base = 10)
     }
     xMin <- min(c(xMin, AllSims[[CurParam]]))
     xMax <- max(c(xMax, AllSims[[CurParam]]))
}
xInitRange <- c(xMin, xMax)
MainEvol <- DispEvol[[CurSpeed]]


FigWidth <- 8
FigHeight <- 6
AllCol <- "lightskyblue"
ExtantCol <- "navyblue"
EvolDispAxisLabel <- expression(paste(Delta, bar(d), sep = ""))
InitDispAxisLabel <- expression(paste("Equilibrium ", italic("log"), "(", italic("d"), ")", sep = ""))
EvolHistBreaks <- seq(-375, 475, by = 5)
InitHistBreaks <- seq(-1.8, 3.4, by = 0.2)
xEvolRange <- c(-100, 100)
Initxticks <- seq(-1, 3, by = 0.2)
AxisSize <- 1.15
OuterMar <- c(0, 1, 0, 0)
InnerMar <- c(5, 5, 4, 3) + 0.1
options(scipen = 1000)
LabSize <- 1.25

pdf(file = "DispComposite.pdf", width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     par(mfrow = c(2,2), mar = InnerMar, oma = OuterMar)
     # DispEvol for parameter 1
     ParamData <- subset(MainEvol, Param == 1)
     ExtantData <- subset(MainEvol, (Param == 1) & (Ext == 1))
     ParamHist <- hist(ParamData$DeltaDisp, plot = FALSE, breaks = EvolHistBreaks)
     plot(ParamHist, col = AllCol, main = "", xlab = EvolDispAxisLabel, ylab = "", 
          las = 1, xlim = xEvolRange, xaxt = "n", yaxt = "n", cex.lab = LabSize)
     axis(1, cex = AxisSize)
     axis(2, cex.axis = AxisSize, las = 1)
     hist(ExtantData$DeltaDisp, col = ExtantCol, breaks = EvolHistBreaks,
          add = TRUE)
     mtext("Frequency", side = 2, line = 4.25, cex = 1.15)
     text(expression(bold("a")), x = 0.9*xEvolRange[1], y = 0.9*max(ParamHist$counts), cex = 1.15)
     legend(legend = c("All sims", "Extant sims"), col = c(AllCol, ExtantCol), pch = 15, cex = 1.15,
            x = "topright", bty = "n")
     # DispEvol for parameter 9
     ParamData <- subset(MainEvol, Param == 9)
     ExtantData <- subset(MainEvol, (Param == 9) & (Ext == 1))
     ParamHist <- hist(ParamData$DeltaDisp, plot = FALSE, breaks = EvolHistBreaks)
     plot(ParamHist, col = AllCol, main = "", xlab = EvolDispAxisLabel, ylab = "", 
          las = 1, xlim = xEvolRange, xaxt = "n", yaxt = "n", cex.lab = LabSize)
     axis(1, cex = AxisSize)
     axis(2, cex.axis = AxisSize, las = 1)
     hist(ExtantData$DeltaDisp, col = ExtantCol, breaks = EvolHistBreaks,
          add = TRUE)
     mtext("Frequency", side = 2, line = 4.25, cex = 1.15)
     text(expression(bold("b")), x = 0.9*xEvolRange[1], y = 0.9*max(ParamHist$counts), cex = 1.15)
     # InitDisp for parameter 1
     AllHist <- hist(AllSims[[1]], plot = FALSE, breaks = InitHistBreaks)
     plot(AllHist, col = AllCol, main = "", xlab = InitDispAxisLabel, ylab = "", cex.axis = AxisSize,
          las = 1, xlim = xInitRange, xaxt = "n", cex.lab = LabSize)
     hist(ExtantSims[[1]], col = ExtantCol, breaks = InitHistBreaks, add = TRUE)
     axis(1, cex.axis = AxisSize)
     axis(1, at = Initxticks, labels = FALSE, tcl = -0.25)
     mtext("Frequency", side = 2, line = 4.25, cex = 1.15)
     text(expression(bold("c")), x = 0.9*xInitRange[1], y = 0.9*max(AllHist$counts), cex = 1.15)
     # InitDisp for parameter 9
     AllHist <- hist(AllSims[[9]], plot = FALSE, breaks = InitHistBreaks)
     plot(AllHist, col = AllCol, main = "", xlab = InitDispAxisLabel, ylab = "", cex.axis = AxisSize,
          las = 1, xlim = xInitRange, xaxt = "n", cex.lab = LabSize)
     hist(ExtantSims[[9]], col = ExtantCol, breaks = InitHistBreaks, add = TRUE)
     axis(1, cex.axis = AxisSize)
     axis(1, at = Initxticks, labels = FALSE, tcl = -0.25)
     mtext("Frequency", side = 2, line = 4.25, cex = 1.15)
     text(expression(bold("d")), x = 0.9*xInitRange[1], y = 0.9*max(AllHist$counts), cex = 1.15)
dev.off()

