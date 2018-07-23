# This script quantifies the probability of individuals with a given dispersal
#    phenotype exactly tracking climate change at three different speeds

log_d <- c(1.2, 1.5, 2)
IdealLoc <- c(1,1,2)
d <- 10 ^ log_d
Nsims <- 100000000
eta <- 50
thetas <- runif(n = Nsims, min = 0, max = 2*pi)
TrackProb <- rep(NA, 3)
# Step through each value, simulate 1 or two rounds of dispersal, and quantify
#    the probability of being in the correct patch to track climate change
for(i in 1:3){
     # Simulate dispersal
     Realized_d <- rexp(n = Nsims, rate = 1 / d[i])
     d_x <- Realized_d * cos(thetas)
     NewPatch <- ifelse(d_x < 0, ceiling((d_x - eta/2) / eta), floor((d_x + eta/2) / eta))
     
     # If we are on the slow speed of climate change, simulate one more round
     #    of dispersal
     if(i == 1){
          NewThetas <- runif(n = Nsims, min = 0, max = 2*pi)
          NewRealized_d <- rexp(n = Nsims, rate = 1 / d[i])
          New_d_x <- NewRealized_d * cos(NewThetas)
          DeltaPatch <- ifelse(New_d_x < 0, ceiling((New_d_x - eta/2) / eta), 
                               floor((New_d_x + eta/2) / eta))
          NewPatch <- NewPatch + DeltaPatch
     }
     
     # Quantify the probability of an individual being in the right patch to track
     #    climate change
     TrackProb[i] <- sum(NewPatch == IdealLoc[i]) / Nsims
}
TrackProb # ~ 9%, 12%, 7%
