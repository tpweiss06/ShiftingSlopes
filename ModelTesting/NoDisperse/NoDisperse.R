# Script to test for normal population growth within a single patch in the 
#    absence of dispersal. 

setwd("~/Desktop/RangeShifts/ShiftingSlopes/ModelTesting/NoDisperse/")
source("~/Desktop/RangeShifts/ShiftingSlopes/SimFunctions.R")

# First set all parameters
BetaInit <- 0
gamma <- 0.025
tau <- 276.25
LocalSel <- 1
omega <- 3
U <- c(0.02, 0)
Vm <- c(0.0004, 0)
names(U) <- names(Vm) <- c("fit", "disp")
nFit <- 5
nDisp <- 5
R0 <- 1
K0 <- 10
width <- 10
kern <- "exp"
EnvGradType <- "K"
monoecious <- TRUE
BurnIn <- 100
BurnOut <- 0
LengthShift <- 0
ClimSpeed <- 0
InitPopSize <- 50
FitInit <- 0.05
FitDiv <- 0.025
DispInit <- -100
DispDiv <- 0
PatchScale <- 10
NumRands <- 1000000
SexRatio <- 0.5


parameters <- list(BetaInit, gamma, tau, LocalSel, omega, U, Vm, nFit, nDisp, R0,
                   K0, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                   LengthShift, ClimSpeed, InitPopSize, FitInit, FitDiv, DispInit,
                   DispDiv, PatchScale, NumRands, SexRatio)
names(parameters) <- c("BetaInit", "gamma", "tau", "LocalSel", "omega", "U", "Vm", 
                       "nFit", "nDisp", "R0", "K0", "width", "kern", "EnvGradType", 
                       "monoecious", "BurnIn", "BurnOut", "LengthShift", "ClimSpeed", 
                       "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                       "PatchScale", "NumRands", "SexRatio")

FullSim(parameters = parameters, parallel = FALSE)


# Graph the resulting abundance patterns through time as a composite graph
TimeVec <- 1:BurnIn
load("Sim1/SummaryStats.Rdata")
PatchAbunds <- PatchAbund[,,57] # rows are generations, columns are patches

pdf(file = "NoDisperseTest.pdf", width = 7.2, height = 5.4)
     par(mfrow = c(2,5), mar = c(0.5,0.5,0.5,0.5), oma = c(c(5,4,1,1)))
     for(p in 1:width){
          plot(TimeVec, PatchAbunds[,p], type = "l", axes = FALSE, xlab = "",
               ylab = "", main = "", lwd = 1.5, ylim = c(0, 140))
          box()
          axis(1, at = seq(0, 100, by = 2), labels = FALSE, tcl = -0.25)
          axis(2, at = seq(0, 140, by = 5), labels = FALSE, tcl = -0.25)
          if(p == 1 | p == 6){
               axis(2, at = seq(0, 140, by = 20), tcl = -0.5, cex.axis = 1.5)
          } else{
               axis(2, at = seq(0, 140, by = 20), labels = FALSE, tcl = -0.5)
          }
          if(p >= 6){
               axis(1, at = seq(0,100,by=10), tcl = -0.5, cex.axis = 1.5)
          } else{
               axis(1, at = seq(0,100,by=10), labels = FALSE, tcl = -0.5)
          }
     }
     mtext("Time", side = 1, outer = TRUE, cex = 1.5, line = 2.5)
     mtext("Abundance", side = 2, outer = TRUE, cex = 1.5, line = 2)
dev.off()




