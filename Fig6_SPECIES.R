# Clear all stored parameters:
rm(list=ls())

##################################################
## Panel A: Extreme underdominance time-series: ##
##################################################

h <- 0.5
s <- 0.1
releaseProportion <- 0.6
mu <- 0.01

qAk <- releaseProportion
pAk <- 1 - qAk

qBk <- 0
pBk <- 1 - qBk

for (k in 1:1000) {
  qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
  qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
  
  pAk[k+1] <- 1 - qAk[k+1]
  pBk[k+1] <- 1 - qBk[k+1]
}

time <- seq(1:40) 
qAk <- qAk[1:40]
pAk <- pAk[1:40]
transgenic_A <- qAk
wildtype_A <- pAk

qBk <- qBk[1:40]
pBk <- pBk[1:40]
transgenic_B <- qBk
wildtype_B <- pBk

trajectory <- cbind(time, qAk, pAk, qBk, pBk,
                    transgenic_A, wildtype_A, transgenic_B, wildtype_B)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

panelA <- ggplot(trajectory_df, aes(x=time, y=trajectory_df, color=Population)) +
  geom_line(aes(y = transgenic_A, col = "A (drive)"), size = 1.2) + 
  geom_line(aes(y = wildtype_A, col = "A (WT)"), size = 1.2) + 
  geom_line(aes(y = transgenic_B, col = "B (drive)"), size = 1.2) + 
  geom_line(aes(y = wildtype_B, col = "B (WT)"), size = 1.2) + 
  labs(x = "Generation", y = "Population frequency") +
  theme(legend.position = c(0.76, 0.5)) + 
  theme(legend.background = element_rect(linetype="solid", color="black")) +
  scale_color_discrete(breaks=c("B (WT)","A (drive)","A (WT)","B (drive)")) + 
  labs(tag = "A")

############################################################################
## Panel B: Extreme underdominance, prevalence in neighboring population: ##
############################################################################

numPoints <- 100
h <- 0.5
s <- 0.1
releaseProportion <- 0.99

muI <- seq(0, 0.2, length = numPoints)
prevPopB <- rep(0, numPoints)

for (i in 1:numPoints) {
  mu <- muI[i]

  qAk <- releaseProportion
  pAk <- 1 - qAk
  
  qBk <- 0
  pBk <- 1 - qBk
  
  for (k in 1:1000) {
    qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
    qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
    
    pAk[k+1] <- 1 - qAk[k+1]
    pBk[k+1] <- 1 - qBk[k+1]
  }
  
  prevPopB[i] <- qBk[1001]
}

prevPopB

trajectory <- cbind(muI, prevPopB)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Migration rate (", italic(mu), ")"))

panelB <- ggplot(trajectory_df, aes(x=muI, y=trajectory_df, col="red")) +
  geom_line(aes(y = prevPopB), size = 1.2, show.legend = FALSE) + 
  labs(x = x_label, y = "Drive frequency (pop B)") +
  coord_cartesian(ylim = c(0, 0.2)) + 
  labs(tag = "B")

##########################################################
## Panel C: Extreme underdominance migration threshold: ##
##########################################################

# Migration thresholds for exterme underdominance, two-population model:

numPoints <- 100
h <- 0.5
sI <- seq(0, 0.8, length = numPoints)
migrationThreshold2Pop <- rep(0, numPoints)
releaseProportion <- 1

for (i in 1:numPoints) {
  s <- sI[i]
  migrationThresholdI <- 0
  
  # Loop to determine region of threshold (0.1):
  for (mu in seq(0, 0.9, by=0.1)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (mu in seq((migrationThresholdI-0.1), migrationThresholdI, by=0.01)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (mu in seq((migrationThresholdI-0.01), migrationThresholdI, by=0.001)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (mu in seq((migrationThresholdI-0.001), migrationThresholdI, by=0.0001)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (mu in seq((migrationThresholdI-0.0001), migrationThresholdI, by=0.00001)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (mu in seq((migrationThresholdI-0.00001), migrationThresholdI, by=0.000001)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  migrationThreshold2Pop[i] <- migrationThresholdI - 0.0000005
}

# Two-population loss thresholds:

migrationThreshold2PopLoss <- rep(0, numPoints)
releaseProportion <- 1

for (i in 1:numPoints) {
  s <- sI[i]
  migrationThresholdI <- 0
  
  # Loop to determine region of threshold (0.1):
  for (mu in seq(0, 0.9, by=0.1)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if ((pBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (mu in seq((migrationThresholdI-0.1), migrationThresholdI, by=0.01)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if ((pBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (mu in seq((migrationThresholdI-0.01), migrationThresholdI, by=0.001)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if ((pBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (mu in seq((migrationThresholdI-0.001), migrationThresholdI, by=0.0001)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if ((pBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (mu in seq((migrationThresholdI-0.0001), migrationThresholdI, by=0.00001)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if ((pBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (mu in seq((migrationThresholdI-0.00001), migrationThresholdI, by=0.000001)) {
    qAk <- releaseProportion
    qBk <- 0
    
    pAk <- 1 - qAk
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s)/((pAk[k]*(1-mu) + pBk[k]*mu)^2 + (qAk[k]*(1-mu) + qBk[k]*mu)^2*(1-s))
      qBk[k+1] <- (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s)/((pBk[k]*(1-mu) + pAk[k]*mu)^2 + (qBk[k]*(1-mu) + qAk[k]*mu)^2*(1-s))
      
      pAk[k+1] <- 1 - qAk[k+1]
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if ((pBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  migrationThreshold2PopLoss[i] <- migrationThresholdI - 0.0000005
}

migrationThreshold2PopLoss

# Migration thresholds for source model:

migrationThresholdSource <- rep(0,numPoints)

for (i in 1:numPoints) {
  s <- sI[i]
  migrationThresholdI <- 0
  
  # First calculate the population A equilibrium:
  releaseProportion <- 0.99
  
  qAk <- releaseProportion
  pAk <- 1 - qAk
  
  for (k in 1:1000) {
    qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
    pAk[k+1] <- 1 - qAk[k+1]
  }
  
  qA <- qAk[1001]
  pA <- pAk[1001]
  
  # Loop to determine region of threshold (0.1):
  for (mu in seq(0, 0.5, by=0.1)) {
    qBk <- 0
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qBk[k+1] <- (qBk[k]*(1-mu) + qA*mu)^2*(1-s)/((pBk[k]*(1-mu) + pA*mu)^2 + (qBk[k]*(1-mu) + qA*mu)^2*(1-s))
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (mu in seq((migrationThresholdI-0.1), migrationThresholdI, by=0.01)) {
    qBk <- 0
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qBk[k+1] <- (qBk[k]*(1-mu) + qA*mu)^2*(1-s)/((pBk[k]*(1-mu) + pA*mu)^2 + (qBk[k]*(1-mu) + qA*mu)^2*(1-s))
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (mu in seq((migrationThresholdI-0.01), migrationThresholdI, by=0.001)) {
    qBk <- 0
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qBk[k+1] <- (qBk[k]*(1-mu) + qA*mu)^2*(1-s)/((pBk[k]*(1-mu) + pA*mu)^2 + (qBk[k]*(1-mu) + qA*mu)^2*(1-s))
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (mu in seq((migrationThresholdI-0.001), migrationThresholdI, by=0.0001)) {
    qBk <- 0
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qBk[k+1] <- (qBk[k]*(1-mu) + qA*mu)^2*(1-s)/((pBk[k]*(1-mu) + pA*mu)^2 + (qBk[k]*(1-mu) + qA*mu)^2*(1-s))
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (mu in seq((migrationThresholdI-0.0001), migrationThresholdI, by=0.00001)) {
    qBk <- 0
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qBk[k+1] <- (qBk[k]*(1-mu) + qA*mu)^2*(1-s)/((pBk[k]*(1-mu) + pA*mu)^2 + (qBk[k]*(1-mu) + qA*mu)^2*(1-s))
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (mu in seq((migrationThresholdI-0.00001), migrationThresholdI, by=0.000001)) {
    qBk <- 0
    pBk <- 1 - qBk
    
    for (k in 1:1000) {
      qBk[k+1] <- (qBk[k]*(1-mu) + qA*mu)^2*(1-s)/((pBk[k]*(1-mu) + pA*mu)^2 + (qBk[k]*(1-mu) + qA*mu)^2*(1-s))
      pBk[k+1] <- 1 - qBk[k+1]
    }
    
    if (qBk[k+1] > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  migrationThresholdSource[i] <- migrationThresholdI - 0.0000005
}

migrationThreshold2Pop[2:numPoints] <- NaN

trajectory <- cbind(sI, migrationThreshold2Pop, migrationThreshold2PopLoss, 
                    migrationThresholdSource)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Fitness cost (", italic("s"), ")"))

panelC <- ggplot(trajectory_df, aes(x=sI, y=trajectory_df, color=Model)) +
  geom_line(aes(y = migrationThreshold2Pop, col = "2 pop"), size = 1.2) + 
  geom_line(aes(y = migrationThreshold2PopLoss, col = "(loss)"), size = 1.2) + 
  geom_line(aes(y = migrationThresholdSource, col = "Source"), size = 1.2) + 
  geom_point(aes(x=sI[1], y=migrationThreshold2Pop[1]), col="#00BA38") +
  coord_cartesian(ylim = c(0, 0.85)) + 
  labs(x = x_label, y = "Migration threshold") + 
  theme(legend.position = c(0.28, 0.66)) + 
  theme(legend.background = element_rect(linetype="solid", color="black")) +
  scale_color_discrete(breaks=c("2 pop","(loss)","Source")) + 
  labs(tag = "C")

#########################################
## Panel D: Translocation time-series: ##
#########################################

h <- 0.5
s <- 0.1
releaseProportion <- 0.6
mu <- 0.01

XXYYAk <- releaseProportion
XxYyAk <- 0
xxyyAk <- 1 - releaseProportion

XXYYBk <- 0
XxYyBk <- 0
xxyyBk <- 1

for (k in 1:1000) {
  XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
  xYA <- 0.25*XxYyAk[k]
  XyA <- 0.25*XxYyAk[k]
  xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
  
  XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
  xYB <- 0.25*XxYyBk[k]
  XyB <- 0.25*XxYyBk[k]
  xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
  
  XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
  XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
  xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
  
  WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
  
  XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
  XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
  xxyyAk[k+1] <- xxyyHatA/WA
  
  XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
  XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
  xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
  
  WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
  
  XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
  XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
  xxyyBk[k+1] <- xxyyHatB/WB
}

time <- seq(1:40) 
XXYYAk <- XXYYAk[1:40]
XxYyAk <- XxYyAk[1:40]
transgenic_A <- XXYYAk + XxYyAk
xxyyAk <- xxyyAk[1:40]
wildtype_A <- xxyyAk

XXYYBk <- XXYYBk[1:40]
XxYyBk <- XxYyBk[1:40]
transgenic_B <- XXYYBk + XxYyBk
xxyyBk <- xxyyBk[1:40]
wildtype_B <- xxyyBk

trajectory <- cbind(time, XXYYAk, XxYyAk, xxyyAk, XXYYBk, XxYyBk, xxyyBk,
                    transgenic_A, wildtype_A, transgenic_B, wildtype_B)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

panelD <- ggplot(trajectory_df, aes(x=time, y=trajectory_df, color=Population)) +
  geom_line(aes(y = transgenic_A, col = "A (drive)"), size = 1.2) + 
  geom_line(aes(y = wildtype_A, col = "A (WT)"), size = 1.2) + 
  geom_line(aes(y = transgenic_B, col = "B (drive)"), size = 1.2) + 
  geom_line(aes(y = wildtype_B, col = "B (WT)"), size = 1.2) + 
  labs(x = "Generation", y = "Population frequency") +
  theme(legend.position = c(0.76, 0.5)) + 
  theme(legend.background = element_rect(linetype="solid", color="black")) +
  scale_color_discrete(breaks=c("B (WT)","A (drive)","A (WT)","B (drive)")) + 
  labs(tag = "D")

####################################################################
## Panel E: Translocations, prevalence in neighboring population: ##
####################################################################

numPoints <- 100
h <- 0.5
s <- 0.1
releaseProportion <- 0.99

muI <- seq(0, 0.2, length = numPoints)
prevPopB <- rep(0, numPoints)

for (i in 1:numPoints) {
  mu <- muI[i]
  
  XXYYAk <- releaseProportion
  XxYyAk <- 0
  xxyyAk <- 1 - releaseProportion
  
  XXYYBk <- 0
  XxYyBk <- 0
  xxyyBk <- 1
  
  for (k in 1:1000) {
    XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
    xYA <- 0.25*XxYyAk[k]
    XyA <- 0.25*XxYyAk[k]
    xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
    
    XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
    xYB <- 0.25*XxYyBk[k]
    XyB <- 0.25*XxYyBk[k]
    xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
    
    XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
    XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
    xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
    
    WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
    
    XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
    XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
    xxyyAk[k+1] <- xxyyHatA/WA
    
    XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
    XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
    xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
    
    WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
    
    XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
    XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
    xxyyBk[k+1] <- xxyyHatB/WB
  }
  
  prevPopB[i] <- 1- xxyyBk[1001]
}

prevPopB

trajectory <- cbind(muI, prevPopB)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Migration rate (", italic(mu), ")"))

panelE <- ggplot(trajectory_df, aes(x=muI, y=trajectory_df, col="red")) +
  geom_line(aes(y = prevPopB), size = 1.2, show.legend = FALSE) + 
  labs(x = x_label, y = "Drive frequency (pop B)") +
  coord_cartesian(ylim = c(0, 0.2)) + 
  labs(tag = "E")

##################################################
## Panel F: Translocations migration threshold: ##
##################################################

# Migration thresholds for translocations, two-population model:

numPoints <- 100
h <- 0.5
sI <- seq(0, 0.8, length = numPoints)
migrationThreshold2Pop <- rep(0, numPoints)
releaseProportion <- 1

for (i in 1:numPoints) {
  s <- sI[i]
  migrationThresholdI <- 0
  
  # Loop to determine region of threshold (0.1):
  for (mu in seq(0, 0.5, by=0.1)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (mu in seq((migrationThresholdI-0.1), migrationThresholdI, by=0.01)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (mu in seq((migrationThresholdI-0.01), migrationThresholdI, by=0.001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (mu in seq((migrationThresholdI-0.001), migrationThresholdI, by=0.0001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (mu in seq((migrationThresholdI-0.0001), migrationThresholdI, by=0.00001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (mu in seq((migrationThresholdI-0.00001), migrationThresholdI, by=0.000001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  migrationThreshold2Pop[i] <- migrationThresholdI - 0.0000005
}

migrationThreshold2Pop

# Two-population model, loss thresholds:

migrationThreshold2PopLoss <- rep(0, numPoints)
releaseProportion <- 1

for (i in 1:numPoints) {
  s <- sI[i]
  migrationThresholdI <- 0
  
  # Loop to determine region of threshold (0.1):
  for (mu in seq(0, 0.5, by=0.1)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (mu in seq((migrationThresholdI-0.1), migrationThresholdI, by=0.01)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (mu in seq((migrationThresholdI-0.01), migrationThresholdI, by=0.001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (mu in seq((migrationThresholdI-0.001), migrationThresholdI, by=0.0001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (mu in seq((migrationThresholdI-0.0001), migrationThresholdI, by=0.00001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (mu in seq((migrationThresholdI-0.00001), migrationThresholdI, by=0.000001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  migrationThreshold2PopLoss[i] <- migrationThresholdI - 0.0000005
}

migrationThreshold2PopLoss

# Migration thresholds for translocations, source model:

migrationThresholdSource <- rep(0,numPoints)

for (i in 1:numPoints) {
  s <- sI[i]
  migrationThresholdI <- 0
  
  # First calculate the population A equilibrium:
  releasePr <- 0.99
  
  XXYYAk <- releasePr
  XxYyAk <- 0
  xxyyAk <- 1-releasePr
  
  for (k in 1:1000) {
    XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
    xYA <- 0.25*XxYyAk[k]
    XyA <- 0.25*XxYyAk[k]
    xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
    
    XXYYHatA <- XYA^2
    XxYyHatA <- 2*(XYA*xyA + xYA*XyA)
    xxyyHatA <- xyA^2
    
    WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
    
    XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
    XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
    xxyyAk[k+1] <- xxyyHatA/WA
  }
  
  XXYYA <- XXYYAk[1001]
  XxYyA <- XxYyAk[1001]
  xxyyA <- xxyyAk[1001]
  
  # Loop to determine region of threshold (0.1):
  for (mu in seq(0, 0.3, by=0.1)) {
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (mu in seq((migrationThresholdI-0.1), migrationThresholdI, by=0.01)) {
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (mu in seq((migrationThresholdI-0.01), migrationThresholdI, by=0.001)) {
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (mu in seq((migrationThresholdI-0.001), migrationThresholdI, by=0.0001)) {
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (mu in seq((migrationThresholdI-0.0001), migrationThresholdI, by=0.00001)) {
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (mu in seq((migrationThresholdI-0.00001), migrationThresholdI, by=0.000001)) {
    XXYYBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.25*XxYyBk[k]
      XyB <- 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  migrationThresholdSource[i] <- migrationThresholdI - 0.0000005
}

migrationThresholdSource[82:100] <- NaN
migrationThreshold2PopLoss[82:100] <- NaN
migrationThreshold2Pop[2:numPoints] <- NaN

trajectory <- cbind(sI, migrationThreshold2Pop, migrationThreshold2PopLoss, 
                    migrationThresholdSource)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Fitness cost (", italic("s"), ")"))

panelF <- ggplot(trajectory_df, aes(x=sI, y=trajectory_df, color=Model)) +
  geom_line(aes(y = migrationThreshold2Pop, col = "2 pop"), size = 1.2) + 
  geom_line(aes(y = migrationThreshold2PopLoss, col = "(loss)"), size = 1.2) + 
  geom_line(aes(y = migrationThresholdSource, col = "Source"), size = 1.2) + 
  geom_point(aes(x=sI[1], y=migrationThreshold2Pop[1]), col="#00BA38") +
  coord_cartesian(ylim = c(0, 0.5)) + 
  labs(x = x_label, y = "Migration threshold") + 
  theme(legend.position = c(0.79, 0.66)) + 
  theme(legend.background = element_rect(linetype="solid", color="black")) +
  scale_color_discrete(breaks=c("2 pop","(loss)","Source")) + 
  labs(tag = "F")

#################################################################
## Panel G: Time-series for 2-locus engineered underdominance: ##
#################################################################

h <- 0.5
s <- 0.1
releaseProportion <- 0.6
mu <- 0.01

XXYYAk <- releaseProportion
XxYYAk <- 0
XXYyAk <- 0
XxYyAk <- 0
xxyyAk <- 1 - releaseProportion

XXYYBk <- 0
XxYYBk <- 0
XXYyBk <- 0
XxYyBk <- 0
xxyyBk <- 1

for (k in 1:1000) {
  XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
  xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
  XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
  xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
  
  XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
  xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
  XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
  xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
  
  XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
  XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
  XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
  XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
  xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
  
  WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
  
  XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
  XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
  XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
  XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
  xxyyAk[k+1] <- xxyyHatA/WA
  
  XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
  XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
  XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
  XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
  xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
  
  WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
  
  XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
  XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
  XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
  XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
  xxyyBk[k+1] <- xxyyHatB/WB
}

time <- seq(1:40) 
XXYYAk <- XXYYAk[1:40]
XxYYAk <- XxYYAk[1:40]
XXYyAk <- XXYyAk[1:40]
XxYyAk <- XxYyAk[1:40]
transgenic_A <- XXYYAk + XxYYAk + XXYyAk + XxYyAk
xxyyAk <- xxyyAk[1:40]
wildtype_A <- xxyyAk

XXYYBk <- XXYYBk[1:40]
XxYYBk <- XxYYBk[1:40]
XXYyBk <- XXYyBk[1:40]
XxYyBk <- XxYyBk[1:40]
transgenic_B <- XXYYBk + XxYYBk + XXYyBk + XxYyBk
xxyyBk <- xxyyBk[1:40]
wildtype_B <- xxyyBk

trajectory <- cbind(time, transgenic_A, wildtype_A, transgenic_B, wildtype_B)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

panelG <- ggplot(trajectory_df, aes(x=time, y=trajectory_df, color=Population)) +
  geom_line(aes(y = transgenic_A, col = "A (drive)"), size = 1.2) + 
  geom_line(aes(y = wildtype_A, col = "A (WT)"), size = 1.2) + 
  geom_line(aes(y = transgenic_B, col = "B (drive)"), size = 1.2) + 
  geom_line(aes(y = wildtype_B, col = "B (WT)"), size = 1.2) + 
  labs(x = "Generation", y = "Population frequency") +
  theme(legend.position = c(0.76, 0.5)) + 
  theme(legend.background = element_rect(linetype="solid", color="black")) +
  scale_color_discrete(breaks=c("B (WT)","A (drive)","A (WT)","B (drive)")) + 
  labs(tag = "G")

#######################################################################################
## Panel H: 2-locus engineered underdominance, prevalence in neighboring population: ##
#######################################################################################

numPoints <- 100
h <- 0.5
s <- 0.1
releaseProportion <- 0.9

muI <- seq(0, 0.2, length = numPoints)
prevPopB <- rep(0, numPoints)

for (i in 1:numPoints) {
  mu <- muI[i]
  
  XXYYAk <- releaseProportion
  XxYYAk <- 0
  XXYyAk <- 0
  XxYyAk <- 0
  xxyyAk <- 1 - releaseProportion
  
  XXYYBk <- 0
  XxYYBk <- 0
  XXYyBk <- 0
  XxYyBk <- 0
  xxyyBk <- 1
  
  for (k in 1:1000) {
    XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
    xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
    XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
    xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
    
    XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
    xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
    XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
    xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
    
    XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
    XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
    XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
    XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
    xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
    
    WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
    
    XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
    XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
    XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
    XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
    xxyyAk[k+1] <- xxyyHatA/WA
    
    XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
    XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
    XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
    XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
    xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
    
    WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
    
    XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
    XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
    XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
    XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
    xxyyBk[k+1] <- xxyyHatB/WB
  }
  
  prevPopB[i] <- 1- xxyyBk[1001]
}

prevPopB

trajectory <- cbind(muI, prevPopB)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

panelH <- ggplot(trajectory_df, aes(x=muI, y=trajectory_df, col="red")) +
  geom_line(aes(y = prevPopB), size = 1.2, show.legend = FALSE) + 
  labs(x = x_label, y = "Drive frequency (pop B)") +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(tag = "H")

#######################################################################
## Panel I: 2-locus engineered underdominance, migration thresholds: ##
#######################################################################

# Migration thresholds for translocations, two-population model:

numPoints <- 100
h <- 0.5
sI <- seq(0, 0.8, length = numPoints)
migrationThreshold2Pop <- rep(0, numPoints)
releaseProportion <- 1

for (i in 1:numPoints) {
  s <- sI[i]
  migrationThresholdI <- 0
  
  # Loop to determine region of threshold (0.1):
  for (mu in seq(0, 0.3, by=0.1)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (mu in seq((migrationThresholdI-0.1), migrationThresholdI, by=0.01)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (mu in seq((migrationThresholdI-0.01), migrationThresholdI, by=0.001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (mu in seq((migrationThresholdI-0.001), migrationThresholdI, by=0.0001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (mu in seq((migrationThresholdI-0.0001), migrationThresholdI, by=0.00001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (mu in seq((migrationThresholdI-0.00001), migrationThresholdI, by=0.000001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  migrationThreshold2Pop[i] <- migrationThresholdI - 0.0000005
}

# Two-population loss thresholds:

migrationThreshold2PopLoss <- rep(0, numPoints)
releaseProportion <- 1

for (i in 1:numPoints) {
  s <- sI[i]
  migrationThresholdI <- 0
  
  # Loop to determine region of threshold (0.1):
  for (mu in seq(0, 0.3, by=0.1)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (mu in seq((migrationThresholdI-0.1), migrationThresholdI, by=0.01)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (mu in seq((migrationThresholdI-0.01), migrationThresholdI, by=0.001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (mu in seq((migrationThresholdI-0.001), migrationThresholdI, by=0.0001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (mu in seq((migrationThresholdI-0.0001), migrationThresholdI, by=0.00001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (mu in seq((migrationThresholdI-0.00001), migrationThresholdI, by=0.000001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatA <- (XYA*(1-mu) + XYB*mu)^2
      XxYYHatA <- 2*(XYA*(1-mu) + XYB*mu)*(xYA*(1-mu) + xYB*mu)
      XXYyHatA <- 2*(XYA*(1-mu) + XYB*mu)*(XyA*(1-mu) + XyB*mu)
      XxYyHatA <- 2*((XYA*(1-mu) + XYB*mu)*(xyA*(1-mu) + xyB*mu)+(xYA*(1-mu) + xYB*mu)*(XyA*(1-mu) + XyB*mu))
      xxyyHatA <- (xyA*(1-mu) + xyB*mu)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
      
      XXYYHatB <- (XYB*(1-mu) + XYA*mu)^2
      XxYYHatB <- 2*(XYB*(1-mu) + XYA*mu)*(xYB*(1-mu) + xYA*mu)
      XXYyHatB <- 2*(XYB*(1-mu) + XYA*mu)*(XyB*(1-mu) + XyA*mu)
      XxYyHatB <- 2*((XYB*(1-mu) + XYA*mu)*(xyB*(1-mu) + xyA*mu)+(xYB*(1-mu) + xYA*mu)*(XyB*(1-mu) + XyA*mu))
      xxyyHatB <- (xyB*(1-mu) + xyA*mu)^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- XXYYHatB*(1-s)/WB
      XxYYBk[k+1] <- XxYYHatB*(1-0.75*s)/WB
      XXYyBk[k+1] <- XXYyHatB*(1-0.75*s)/WB
      XxYyBk[k+1] <- XxYyHatB*(1-0.5*s)/WB
      xxyyBk[k+1] <- xxyyHatB/WB
    }
    
    if ((xxyyBk[k+1] == 1) && (mu > 0)) {
      migrationThresholdI <- mu
      break
    }
  }
  
  migrationThreshold2PopLoss[i] <- migrationThresholdI - 0.0000005
}

migrationThreshold2PopLoss

# Migration thresholds for source model:

migrationThresholdSource <- rep(0,numPoints)

for (i in 1:numPoints) {
  s <- sI[i]
  migrationThresholdI <- 0
  
  # First calculate the population A equilibrium:
  releasePr <- 0.99
  
  XXYYAk <- releasePr
  XxYYAk <- 0
  XXYyAk <- 0
  XxYyAk <- 0
  xxyyAk <- 1-releasePr
  
  for (k in 1:1000) {
    XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
    xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
    XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
    xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
    
    XXYYHatA <- XYA^2
    XxYYHatA <- 2*XYA*xYA
    XXYyHatA <- 2*XYA*XyA
    XxYyHatA <- 2*(XYA*xyA + xYA*XyA)
    xxyyHatA <- xyA^2
    
    WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
    
    XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
    XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
    XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
    XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
    xxyyAk[k+1] <- xxyyHatA/WA
  }
  
  XXYYA <- XXYYAk[1001]
  XxYYA <- XxYYAk[1001]
  XXYyA <- XXYyAk[1001]
  XxYyA <- XxYyAk[1001]
  xxyyA <- xxyyAk[1001]
  
  # Loop to determine region of threshold (0.1):
  for (mu in seq(0, 0.3, by=0.1)) {
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYYHatB <- 2*XYB*xYB
      XXYyHatB <- 2*XYB*XyB
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYYBk[k+1] <- (XxYYHatB*(1-0.75*s) + XxYYA*mu)/(WB + mu)
      XXYyBk[k+1] <- (XXYyHatB*(1-0.75*s) + XXYyA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (mu in seq((migrationThresholdI-0.1), migrationThresholdI, by=0.01)) {
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYYHatB <- 2*XYB*xYB
      XXYyHatB <- 2*XYB*XyB
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYYBk[k+1] <- (XxYYHatB*(1-0.75*s) + XxYYA*mu)/(WB + mu)
      XXYyBk[k+1] <- (XXYyHatB*(1-0.75*s) + XXYyA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (mu in seq((migrationThresholdI-0.01), migrationThresholdI, by=0.001)) {
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYYHatB <- 2*XYB*xYB
      XXYyHatB <- 2*XYB*XyB
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYYBk[k+1] <- (XxYYHatB*(1-0.75*s) + XxYYA*mu)/(WB + mu)
      XXYyBk[k+1] <- (XXYyHatB*(1-0.75*s) + XXYyA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (mu in seq((migrationThresholdI-0.001), migrationThresholdI, by=0.0001)) {
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYYHatB <- 2*XYB*xYB
      XXYyHatB <- 2*XYB*XyB
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYYBk[k+1] <- (XxYYHatB*(1-0.75*s) + XxYYA*mu)/(WB + mu)
      XXYyBk[k+1] <- (XXYyHatB*(1-0.75*s) + XXYyA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (mu in seq((migrationThresholdI-0.0001), migrationThresholdI, by=0.00001)) {
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYYHatB <- 2*XYB*xYB
      XXYyHatB <- 2*XYB*XyB
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYYBk[k+1] <- (XxYYHatB*(1-0.75*s) + XxYYA*mu)/(WB + mu)
      XXYyBk[k+1] <- (XXYyHatB*(1-0.75*s) + XXYyA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (mu in seq((migrationThresholdI-0.00001), migrationThresholdI, by=0.000001)) {
    XXYYBk <- 0
    XxYYBk <- 0
    XXYyBk <- 0
    XxYyBk <- 0
    xxyyBk <- 1
    
    for (k in 1:1000) {
      XYB <- XXYYBk[k] + 0.5*XxYYBk[k] + 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xYB <- 0.5*XxYYBk[k] + 0.25*XxYyBk[k]
      XyB <- 0.5*XXYyBk[k] + 0.25*XxYyBk[k]
      xyB <- xxyyBk[k] + 0.25*XxYyBk[k]
      
      XXYYHatB <- XYB^2
      XxYYHatB <- 2*XYB*xYB
      XXYyHatB <- 2*XYB*XyB
      XxYyHatB <- 2*(XYB*xyB + xYB*XyB)
      xxyyHatB <- xyB^2
      
      WB <- XXYYHatB*(1-s) + (XxYYHatB+XXYyHatB)*(1-0.75*s) + XxYyHatB*(1-0.5*s) + xxyyHatB
      
      XXYYBk[k+1] <- (XXYYHatB*(1-s) + XXYYA*mu)/(WB + mu)
      XxYYBk[k+1] <- (XxYYHatB*(1-0.75*s) + XxYYA*mu)/(WB + mu)
      XXYyBk[k+1] <- (XXYyHatB*(1-0.75*s) + XXYyA*mu)/(WB + mu)
      XxYyBk[k+1] <- (XxYyHatB*(1-0.5*s) + XxYyA*mu)/(WB + mu)
      xxyyBk[k+1] <- (xxyyHatB + xxyyA*mu)/(WB + mu)
    }
    
    if ((1-xxyyBk[k+1]) > 0.5) {
      migrationThresholdI <- mu
      break
    }
  }
  
  migrationThresholdSource[i] <- migrationThresholdI - 0.0000005
}

migrationThresholdSource[92:100] <- NaN
migrationThreshold2PopLoss[92:100] <- NaN
migrationThreshold2PopLoss[1:45] <- NaN
migrationThreshold2Pop[52:100] <- NaN

trajectory <- cbind(sI, migrationThreshold2Pop, migrationThreshold2PopLoss, 
                    migrationThresholdSource)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Fitness cost (", italic("s"), ")"))

panelI <- ggplot(trajectory_df, aes(x=sI, y=trajectory_df, color=Model)) +
  geom_line(aes(y = migrationThreshold2Pop, col = "2 pop"), size = 1.2) + 
  geom_line(aes(y = migrationThreshold2PopLoss, col = "(loss)"), size = 1.2) + 
  geom_line(aes(y = migrationThresholdSource, col = "Source"), size = 1.2) + 
  coord_cartesian(ylim = c(0, 0.5)) + 
  labs(x = x_label, y = "Migration threshold") + 
  theme(legend.position = c(0.79, 0.66)) + 
  theme(legend.background = element_rect(linetype="solid", color="black")) +
  scale_color_discrete(breaks=c("2 pop","(loss)","Source")) + 
  labs(tag = "I")

#######################################
## Combine all panels into one plot: ##
#######################################

library(gridExtra)

tiff("test.tiff", units="in", width=8.7, height=8.7*(680/786), res=300)

grid.arrange(panelA, panelB, panelC, panelD, panelE, panelF, 
             panelG, panelH, panelI, ncol=3, nrow=3) # Combine the plots

dev.off()
