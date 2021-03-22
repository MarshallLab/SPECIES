# Clear all stored parameters:
rm(list=ls())

##################################################
## Panel A: Extreme underdominance time-series: ##
##################################################

h <- 0.5
s <- 0.1

# 45% release frequency:
releaseProportion <- 0.45
qAk <- releaseProportion
pAk <- 1 - qAk
for (k in 1:50) {
  qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
  pAk[k+1] <- 1 - qAk[k+1]
}
time <- seq(1:30) 
qAk <- qAk[1:30]
pAk <- pAk[1:30]
transgenic_45 <- qAk

# 50% release frequency:
releaseProportion <- 0.50
qAk <- releaseProportion
pAk <- 1 - qAk
for (k in 1:50) {
  qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
  pAk[k+1] <- 1 - qAk[k+1]
}
time <- seq(1:30) 
qAk <- qAk[1:30]
pAk <- pAk[1:30]
transgenic_50 <- qAk

# 55% release frequency:
releaseProportion <- 0.55
qAk <- releaseProportion
pAk <- 1 - qAk
for (k in 1:50) {
  qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
  pAk[k+1] <- 1 - qAk[k+1]
}
time <- seq(1:30) 
qAk <- qAk[1:30]
pAk <- pAk[1:30]
transgenic_55 <- qAk

# 60% release frequency:
releaseProportion <- 0.60
qAk <- releaseProportion
pAk <- 1 - qAk
for (k in 1:50) {
  qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
  pAk[k+1] <- 1 - qAk[k+1]
}
time <- seq(1:30) 
qAk <- qAk[1:30]
pAk <- pAk[1:30]
transgenic_60 <- qAk

trajectory <- cbind(time, transgenic_45, transgenic_50, transgenic_55, transgenic_60)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

panelA <- ggplot(trajectory_df, aes(x=time, y=trajectory_df, color=Release)) +
  geom_line(aes(y = transgenic_45, col = "45%"), size = 1.2) + 
  geom_line(aes(y = transgenic_50, col = "50%"), size = 1.2) + 
  geom_line(aes(y = transgenic_55, col = "55%"), size = 1.2) + 
  geom_line(aes(y = transgenic_60, col = "60%"), size = 1.2) + 
  labs(x = "Generation", y = "Population frequency") +
  theme(legend.position = c(0.81, 0.51)) + 
  theme(legend.background = element_rect(linetype="solid", color="black")) +
  scale_color_discrete(breaks=c("60%","55%","50%","45%")) + 
  labs(tag = "A")

########################################################
## Panel B: Extreme underdominance release threshold: ##
########################################################

numPoints <- 100
h <- 0.5
sI <- seq(0, 0.8, length = numPoints)
releaseThreshold <- rep(0, numPoints)

for (i in 1:numPoints) {
  s <- sI[i]
  releaseThresholdI <- 0
  
  # Loop to determine region of threshold (0.1):
  for (releaseProportion in seq(0, 1, by=0.1)) {
    qAk <- releaseProportion
    pAk <- 1 - qAk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
      pAk[k+1] <- 1 - qAk[k+1]
    }
    
    if (qAk[k+1] > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (releaseProportion in seq((releaseThresholdI-0.1), releaseThresholdI, by=0.01)) {
    qAk <- releaseProportion
    pAk <- 1 - qAk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
      pAk[k+1] <- 1 - qAk[k+1]
    }
    
    if (qAk[k+1] > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (releaseProportion in seq((releaseThresholdI-0.01), releaseThresholdI, by=0.001)) {
    qAk <- releaseProportion
    pAk <- 1 - qAk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
      pAk[k+1] <- 1 - qAk[k+1]
    }
    
    if (qAk[k+1] > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (releaseProportion in seq((releaseThresholdI-0.001), releaseThresholdI, by=0.0001)) {
    qAk <- releaseProportion
    pAk <- 1 - qAk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
      pAk[k+1] <- 1 - qAk[k+1]
    }
    
    if (qAk[k+1] > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (releaseProportion in seq((releaseThresholdI-0.0001), releaseThresholdI, by=0.00001)) {
    qAk <- releaseProportion
    pAk <- 1 - qAk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
      pAk[k+1] <- 1 - qAk[k+1]
    }
    
    if (qAk[k+1] > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (releaseProportion in seq((releaseThresholdI-0.00001), releaseThresholdI, by=0.000001)) {
    qAk <- releaseProportion
    pAk <- 1 - qAk
    
    for (k in 1:1000) {
      qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
      pAk[k+1] <- 1 - qAk[k+1]
    }
    
    if (qAk[k+1] > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  releaseThreshold[i] <- releaseThresholdI - 0.0000005
}

releaseThreshold
sI

trajectory <- cbind(sI, releaseThreshold)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Fitness cost (", italic("s"), ")"))

panelB <- ggplot(trajectory_df, aes(x=sI, y=trajectory_df, col="red")) +
  geom_line(aes(y = releaseThreshold), size = 1.2, show.legend = FALSE) + 
  labs(x = x_label, y = "Release threshold") + 
  coord_cartesian(ylim = c(0.25, 1)) + 
  labs(tag = "B")

############################################################
## Panel C: Time taken to reach 99% transgenic frequency: ##
############################################################

numPoints <- 100
sI <- seq(0, 0.2, length = numPoints)
timeTo99 <- rep(0, numPoints)
releaseProportion <- 0.65

for (i in 1:numPoints) {
  s <- sI[i]
  timeTo99I <- 0
  
  qAk <- releaseProportion
  pAk <- 1 - qAk
  
  for (k in 1:1000) {
    qAk[k+1] <- (qAk[k])^2*(1-s)/((pAk[k])^2 + (qAk[k])^2*(1-s))
    pAk[k+1] <- 1 - qAk[k+1]
    
    if (pAk[k+1] < 0.01) {
      timeTo99I <- k
      break
    }
  }
  
  timeTo99[i] <- timeTo99I
}

trajectory <- cbind(sI, timeTo99)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Fitness cost (", italic("s"), ")"))

panelC <- ggplot(trajectory_df, aes(x=sI, y=trajectory_df, col="red")) +
  geom_line(aes(y = timeTo99), size = 1.2, show.legend = FALSE) + 
  labs(x = x_label, y = "Generations to 99%") +
  coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 15)) + 
  labs(tag = "C")

#########################################
## Panel D: Translocation time-series: ##
#########################################

h <- 0.5
s <- 0.1

# 50% release frequency:
releaseProportion <- 0.50
XXYYAk <- releaseProportion
XxYyAk <- 0
xxyyAk <- 1 - releaseProportion

for (k in 1:50) {
  XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
  xYA <- 0.25*XxYyAk[k]
  XyA <- 0.25*XxYyAk[k]
  xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
  
  XXYYHatA <- (XYA)^2
  XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
  xxyyHatA <- (xyA)^2
  
  WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
  
  XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
  XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
  xxyyAk[k+1] <- xxyyHatA/WA
}
time <- seq(1:30) 
xxyyA <- xxyyAk[1:30]
transgenic_50 <- 1 - xxyyA

# 55% release frequency:
releaseProportion <- 0.55
XXYYAk <- releaseProportion
XxYyAk <- 0
xxyyAk <- 1 - releaseProportion

for (k in 1:50) {
  XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
  xYA <- 0.25*XxYyAk[k]
  XyA <- 0.25*XxYyAk[k]
  xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
  
  XXYYHatA <- (XYA)^2
  XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
  xxyyHatA <- (xyA)^2
  
  WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
  
  XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
  XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
  xxyyAk[k+1] <- xxyyHatA/WA
}
time <- seq(1:30) 
xxyyA <- xxyyAk[1:30]
transgenic_55 <- 1 - xxyyA

# 60% release frequency:
releaseProportion <- 0.60
XXYYAk <- releaseProportion
XxYyAk <- 0
xxyyAk <- 1 - releaseProportion

for (k in 1:50) {
  XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
  xYA <- 0.25*XxYyAk[k]
  XyA <- 0.25*XxYyAk[k]
  xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
  
  XXYYHatA <- (XYA)^2
  XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
  xxyyHatA <- (xyA)^2
  
  WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
  
  XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
  XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
  xxyyAk[k+1] <- xxyyHatA/WA
}
time <- seq(1:30) 
xxyyA <- xxyyAk[1:30]
transgenic_60 <- 1 - xxyyA

# 65% release frequency:
releaseProportion <- 0.65
XXYYAk <- releaseProportion
XxYyAk <- 0
xxyyAk <- 1 - releaseProportion

for (k in 1:50) {
  XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
  xYA <- 0.25*XxYyAk[k]
  XyA <- 0.25*XxYyAk[k]
  xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
  
  XXYYHatA <- (XYA)^2
  XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
  xxyyHatA <- (xyA)^2
  
  WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
  
  XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
  XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
  xxyyAk[k+1] <- xxyyHatA/WA
}
time <- seq(1:30) 
xxyyA <- xxyyAk[1:30]
transgenic_65 <- 1 - xxyyA

trajectory <- cbind(time, transgenic_50, transgenic_55, transgenic_60, transgenic_65)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

panelD <- ggplot(trajectory_df, aes(x=time, y=trajectory_df, color=Release)) +
  geom_line(aes(y = transgenic_50, col = "50%"), size = 1.2) + 
  geom_line(aes(y = transgenic_55, col = "55%"), size = 1.2) + 
  geom_line(aes(y = transgenic_60, col = "60%"), size = 1.2) + 
  geom_line(aes(y = transgenic_65, col = "65%"), size = 1.2) + 
  labs(x = "Generation", y = "Population frequency") +
  theme(legend.position = c(0.81, 0.51)) + 
  theme(legend.background = element_rect(linetype="solid", color="black")) +
  scale_color_discrete(breaks=c("65%","60%","55%","50%")) + 
  labs(tag = "D")

#####################################################
## Panel E: Release thresholds for translocations: ##
#####################################################

numPoints <- 100
h <- 0.5
sI <- seq(0, 0.8, length = numPoints)
releaseThreshold <- rep(0, numPoints)

for (i in 1:numPoints) {
  s <- sI[i]
  releaseThresholdI <- 0
  
  # Loop to determine region of threshold (0.1):
  for (releaseProportion in seq(0, 1, by=0.1)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion

    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]

      XXYYHatA <- (XYA)^2
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (releaseProportion in seq((releaseThresholdI-0.1), releaseThresholdI, by=0.01)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (releaseProportion in seq((releaseThresholdI-0.01), releaseThresholdI, by=0.001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (releaseProportion in seq((releaseThresholdI-0.001), releaseThresholdI, by=0.0001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (releaseProportion in seq((releaseThresholdI-0.0001), releaseThresholdI, by=0.00001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (releaseProportion in seq((releaseThresholdI-0.00001), releaseThresholdI, by=0.000001)) {
    XXYYAk <- releaseProportion
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.25*XxYyAk[k]
      XyA <- 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  releaseThreshold[i] <- releaseThresholdI - 0.0000005
}

releaseThreshold
sI

trajectory <- cbind(sI, releaseThreshold)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Fitness cost (", italic("s"), ")"))

panelE <- ggplot(trajectory_df, aes(x=sI, y=trajectory_df, col="red")) +
  geom_line(aes(y = releaseThreshold), size = 1.2, show.legend = FALSE) + 
  labs(x = x_label, y = "Release threshold") + 
  coord_cartesian(ylim = c(0.25, 1)) + 
  labs(tag = "E")

############################################################
## Panel F: Time taken to reach 99% transgenic frequency: ##
############################################################

numPoints <- 100
h <- 0.5
sI <- seq(0, 0.20, length = numPoints)
timeTo99 <- rep(0, numPoints)
releaseProportion <- 0.65

for (i in 1:numPoints) {
  s <- sI[i]
  timeTo99I <- 0
  
  XXYYAk <- releaseProportion
  XxYyAk <- 0
  xxyyAk <- 1 - releaseProportion
  
  for (k in 1:1000) {
    XYA <- XXYYAk[k] + 0.25*XxYyAk[k]
    xYA <- 0.25*XxYyAk[k]
    XyA <- 0.25*XxYyAk[k]
    xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
    
    XXYYHatA <- (XYA)^2
    XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
    xxyyHatA <- (xyA)^2
    
    WA <- XXYYHatA*(1-s) + XxYyHatA*(1-0.5*s) + xxyyHatA
    
    XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
    XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
    xxyyAk[k+1] <- xxyyHatA/WA
    
    if (xxyyAk[k+1] < 0.01) {
      timeTo99I <- k
      break
    }
  }
  
  timeTo99[i] <- timeTo99I
}

trajectory <- cbind(sI, timeTo99)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Fitness cost (", italic("s"), ")"))

panelF <- ggplot(trajectory_df, aes(x=sI, y=trajectory_df, col="red")) +
  geom_line(aes(y = timeTo99), size = 1.2, show.legend = FALSE) + 
  labs(x = x_label, y = "Generations to 99%") +
  coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 15)) + 
  labs(tag = "F")

#################################################################
## Panel G: Time-series for 2-locus engineered underdominance: ##
#################################################################

h <- 0.5
s <- 0.1

# 25% release frequency:
releaseProportion <- 0.25
XXYYAk <- releaseProportion
XxYYAk <- 0
XXYyAk <- 0
XxYyAk <- 0
xxyyAk <- 1 - releaseProportion

for (k in 1:1000) {
  XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
  xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
  XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
  xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
  
  XXYYHatA <- (XYA)^2
  XxYYHatA <- 2*(XYA)*(xYA)
  XXYyHatA <- 2*(XYA)*(XyA)
  XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
  xxyyHatA <- (xyA)^2
  
  WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
  
  XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
  XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
  XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
  XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
  xxyyAk[k+1] <- xxyyHatA/WA
}
time <- seq(1:30) 
xxyyA <- xxyyAk[1:30]
transgenic_25 <- 1 - xxyyA

# 30% release frequency:
releaseProportion <- 0.30
XXYYAk <- releaseProportion
XxYYAk <- 0
XXYyAk <- 0
XxYyAk <- 0
xxyyAk <- 1 - releaseProportion

for (k in 1:1000) {
  XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
  xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
  XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
  xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
  
  XXYYHatA <- (XYA)^2
  XxYYHatA <- 2*(XYA)*(xYA)
  XXYyHatA <- 2*(XYA)*(XyA)
  XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
  xxyyHatA <- (xyA)^2
  
  WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
  
  XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
  XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
  XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
  XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
  xxyyAk[k+1] <- xxyyHatA/WA
}
time <- seq(1:30) 
xxyyA <- xxyyAk[1:30]
transgenic_30 <- 1 - xxyyA

# 35% release frequency:
releaseProportion <- 0.35
XXYYAk <- releaseProportion
XxYYAk <- 0
XXYyAk <- 0
XxYyAk <- 0
xxyyAk <- 1 - releaseProportion

for (k in 1:1000) {
  XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
  xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
  XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
  xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
  
  XXYYHatA <- (XYA)^2
  XxYYHatA <- 2*(XYA)*(xYA)
  XXYyHatA <- 2*(XYA)*(XyA)
  XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
  xxyyHatA <- (xyA)^2
  
  WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
  
  XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
  XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
  XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
  XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
  xxyyAk[k+1] <- xxyyHatA/WA
}
time <- seq(1:30) 
xxyyA <- xxyyAk[1:30]
transgenic_35 <- 1 - xxyyA

# 40% release frequency:
releaseProportion <- 0.40
XXYYAk <- releaseProportion
XxYYAk <- 0
XXYyAk <- 0
XxYyAk <- 0
xxyyAk <- 1 - releaseProportion

for (k in 1:1000) {
  XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
  xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
  XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
  xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
  
  XXYYHatA <- (XYA)^2
  XxYYHatA <- 2*(XYA)*(xYA)
  XXYyHatA <- 2*(XYA)*(XyA)
  XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
  xxyyHatA <- (xyA)^2
  
  WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
  
  XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
  XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
  XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
  XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
  xxyyAk[k+1] <- xxyyHatA/WA
}
time <- seq(1:30) 
xxyyA <- xxyyAk[1:30]
transgenic_40 <- 1 - xxyyA

trajectory <- cbind(time, transgenic_40, transgenic_35, transgenic_30, transgenic_25)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

panelG <- ggplot(trajectory_df, aes(x=time, y=trajectory_df, color=Release)) +
  geom_line(aes(y = transgenic_40, col = "40%"), size = 1.2) + 
  geom_line(aes(y = transgenic_35, col = "35%"), size = 1.2) + 
  geom_line(aes(y = transgenic_30, col = "30%"), size = 1.2) + 
  geom_line(aes(y = transgenic_25, col = "25%"), size = 1.2) + 
  labs(x = "Generation", y = "Population frequency") +
  theme(legend.position = c(0.81, 0.51)) + 
  theme(legend.background = element_rect(linetype="solid", color="black")) +
  scale_color_discrete(breaks=c("40%","35%","30%","25%")) + 
  labs(tag = "G")

########################################################################
## Panel H: Release thresholds for 2-locus engineered underdominance: ##
########################################################################

numPoints <- 100
h <- 0.5
sI <- seq(0, 0.8, length = numPoints)
releaseThreshold <- rep(0, numPoints)

for (i in 1:numPoints) {
  s <- sI[i]
  releaseThresholdI <- 0
  
  # Loop to determine region of threshold (0.1):
  for (releaseProportion in seq(0, 1, by=0.1)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYYHatA <- 2*(XYA)*(xYA)
      XXYyHatA <- 2*(XYA)*(XyA)
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.01):
  for (releaseProportion in seq((releaseThresholdI-0.1), releaseThresholdI, by=0.01)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYYHatA <- 2*(XYA)*(xYA)
      XXYyHatA <- 2*(XYA)*(XyA)
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.001):
  for (releaseProportion in seq((releaseThresholdI-0.01), releaseThresholdI, by=0.001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYYHatA <- 2*(XYA)*(xYA)
      XXYyHatA <- 2*(XYA)*(XyA)
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.0001):
  for (mu in seq((releaseThresholdI-0.001), releaseThresholdI, by=0.0001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYYHatA <- 2*(XYA)*(xYA)
      XXYyHatA <- 2*(XYA)*(XyA)
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.00001):
  for (releaseProportion in seq((releaseThresholdI-0.0001), releaseThresholdI, by=0.00001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYYHatA <- 2*(XYA)*(xYA)
      XXYyHatA <- 2*(XYA)*(XyA)
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  # Loop to narrow range of threshold (0.000001):
  for (releaseProportion in seq((releaseThresholdI-0.00001), releaseThresholdI, by=0.000001)) {
    XXYYAk <- releaseProportion
    XxYYAk <- 0
    XXYyAk <- 0
    XxYyAk <- 0
    xxyyAk <- 1 - releaseProportion
    
    for (k in 1:1000) {
      XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
      XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
      xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
      
      XXYYHatA <- (XYA)^2
      XxYYHatA <- 2*(XYA)*(xYA)
      XXYyHatA <- 2*(XYA)*(XyA)
      XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
      xxyyHatA <- (xyA)^2
      
      WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
      
      XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
      XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
      XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
      XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
      xxyyAk[k+1] <- xxyyHatA/WA
    }
    
    if ((1-xxyyAk[k+1]) > 0.9) {
      releaseThresholdI <- releaseProportion
      break
    }
  }
  
  releaseThreshold[i] <- releaseThresholdI - 0.0000005
}

releaseThreshold
sI

trajectory <- cbind(sI, releaseThreshold)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Fitness cost (", italic("s"), ")"))

panelH <- ggplot(trajectory_df, aes(x=sI, y=trajectory_df, col="red")) +
  geom_line(aes(y = releaseThreshold), size = 1.2, show.legend = FALSE) + 
  labs(x = x_label, y = "Release threshold") + 
  coord_cartesian(ylim = c(0.25, 1)) + 
  labs(tag = "H")

############################################################
## Panel I: Time taken to reach 99% transgenic frequency: ##
############################################################

numPoints <- 100
sI <- seq(0, 0.20, length = numPoints)
timeTo99 <- rep(0, numPoints)
releaseProportion <- 0.65

for (i in 1:numPoints) {
  s <- sI[i]
  timeTo99I <- 0
  
  XXYYAk <- releaseProportion
  XxYYAk <- 0
  XXYyAk <- 0
  XxYyAk <- 0
  xxyyAk <- 1 - releaseProportion
  
  for (k in 1:1000) {
    XYA <- XXYYAk[k] + 0.5*XxYYAk[k] + 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
    xYA <- 0.5*XxYYAk[k] + 0.25*XxYyAk[k]
    XyA <- 0.5*XXYyAk[k] + 0.25*XxYyAk[k]
    xyA <- xxyyAk[k] + 0.25*XxYyAk[k]
    
    XXYYHatA <- (XYA)^2
    XxYYHatA <- 2*(XYA)*(xYA)
    XXYyHatA <- 2*(XYA)*(XyA)
    XxYyHatA <- 2*((XYA)*(xyA)+(xYA)*(XyA))
    xxyyHatA <- (xyA)^2
    
    WA <- XXYYHatA*(1-s) + (XxYYHatA+XXYyHatA)*(1-0.75*s) + XxYyHatA*(1-0.5*s) + xxyyHatA
    
    XXYYAk[k+1] <- XXYYHatA*(1-s)/WA
    XxYYAk[k+1] <- XxYYHatA*(1-0.75*s)/WA
    XXYyAk[k+1] <- XXYyHatA*(1-0.75*s)/WA
    XxYyAk[k+1] <- XxYyHatA*(1-0.5*s)/WA
    xxyyAk[k+1] <- xxyyHatA/WA
    
    if (xxyyAk[k+1] < 0.01) {
      timeTo99I <- k
      break
    }
  }
  
  timeTo99[i] <- timeTo99I
}

timeTo99

trajectory <- cbind(sI, timeTo99)
trajectory_df <- data.frame(trajectory)
trajectory
trajectory_df

library(ggplot2)

x_label <- expression(paste("Fitness cost (", italic("s"), ")"))

panelI <- ggplot(trajectory_df, aes(x=sI, y=trajectory_df, col="red")) +
  geom_line(aes(y = timeTo99), size = 1.2, show.legend = FALSE) + 
  labs(x = x_label, y = "Generations to 99%") +
  coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 15)) + 
  labs(tag = "I")

#######################################
## Combine all panels into one plot: ##
#######################################

library(gridExtra)

tiff("test.tiff", units="in", width=8.7, height=8.7*(680/786), res=300)

grid.arrange(panelA, panelB, panelC, panelD, panelE, panelF, 
             panelG, panelH, panelI, ncol=3, nrow=3) # Combine the plots

dev.off()
