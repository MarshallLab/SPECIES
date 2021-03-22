#########################################################################
## SPECIES model fit (construct A1)                                    ##
## Code: John Marshall (john.marshall@berkeley.edu)                    ##
######################################################################### 

# Clear all stored parameters:
rm(list=ls())

# Load required packages:
library(car)
library(truncnorm)
library(ggplot2)
library(coda)
library(grid)

#########################################################################
## Load the data:                                                      ##
#########################################################################

## Read in the data from drive experiments:
setwd("Insert directory where your data is here...")

Data <- read.csv("UnderdominanceData.csv")
Data
head(Data)

gens <- seq(0,11)

# Data for construct 1045A:
X1045A_Rep1_GFPp <- Data$X1045A_Rep1_GFPp
X1045A_Rep1_GFPm <- Data$X1045A_Rep1_GFPm
X1045A_Rep2_GFPp <- Data$X1045A_Rep2_GFPp
X1045A_Rep2_GFPm <- Data$X1045A_Rep2_GFPm
X1045A_Rep3_GFPp <- Data$X1045A_Rep3_GFPp
X1045A_Rep3_GFPm <- Data$X1045A_Rep3_GFPm 
X1045A_Rep4_GFPp <- Data$X1045A_Rep4_GFPp 
X1045A_Rep4_GFPm <- Data$X1045A_Rep4_GFPm
X1045A_Rep5_GFPp <- Data$X1045A_Rep5_GFPp 
X1045A_Rep5_GFPm <- Data$X1045A_Rep5_GFPm
X1045A_Rep6_GFPp <- Data$X1045A_Rep6_GFPp 
X1045A_Rep6_GFPm <- Data$X1045A_Rep6_GFPm 
X1045A_Rep7_GFPp <- Data$X1045A_Rep7_GFPp 
X1045A_Rep7_GFPm <- Data$X1045A_Rep7_GFPm 
X1045A_Rep8_GFPp <- Data$X1045A_Rep8_GFPp
X1045A_Rep8_GFPm <- Data$X1045A_Rep8_GFPm 
X1045A_Rep9_GFPp <- Data$X1045A_Rep9_GFPp 
X1045A_Rep9_GFPm <- Data$X1045A_Rep9_GFPm 

# Data for construct A1:
A1_Rep1_GFPp <- Data$A1_Rep1_GFPp 
A1_Rep1_GFPm <- Data$A1_Rep1_GFPm 
A1_Rep2_GFPp <- Data$A1_Rep2_GFPp 
A1_Rep2_GFPm <- Data$A1_Rep2_GFPm 
A1_Rep3_GFPp <- Data$A1_Rep3_GFPp 
A1_Rep3_GFPm <- Data$A1_Rep3_GFPm 
A1_Rep4_GFPp <- Data$A1_Rep4_GFPp 
A1_Rep4_GFPm <- Data$A1_Rep4_GFPm
A1_Rep5_GFPp <- Data$A1_Rep5_GFPp 
A1_Rep5_GFPm <- Data$A1_Rep5_GFPm 

# Totals for construct 1045A:
Total_1045A_Rep1 <- (X1045A_Rep1_GFPp + X1045A_Rep1_GFPm)
Total_1045A_Rep2 <- (X1045A_Rep2_GFPp + X1045A_Rep2_GFPm)
Total_1045A_Rep3 <- (X1045A_Rep3_GFPp + X1045A_Rep3_GFPm)
Total_1045A_Rep4 <- (X1045A_Rep4_GFPp + X1045A_Rep4_GFPm)
Total_1045A_Rep5 <- (X1045A_Rep5_GFPp + X1045A_Rep5_GFPm)
Total_1045A_Rep6 <- (X1045A_Rep6_GFPp + X1045A_Rep6_GFPm)
Total_1045A_Rep7 <- (X1045A_Rep7_GFPp + X1045A_Rep7_GFPm)
Total_1045A_Rep8 <- (X1045A_Rep8_GFPp + X1045A_Rep8_GFPm)
Total_1045A_Rep9 <- (X1045A_Rep9_GFPp + X1045A_Rep9_GFPm)

# Totals for construct A1:
Total_A1_Rep1 <- (A1_Rep1_GFPp + A1_Rep1_GFPm)
Total_A1_Rep2 <- (A1_Rep2_GFPp + A1_Rep2_GFPm)
Total_A1_Rep3 <- (A1_Rep3_GFPp + A1_Rep3_GFPm)
Total_A1_Rep4 <- (A1_Rep4_GFPp + A1_Rep4_GFPm)
Total_A1_Rep5 <- (A1_Rep5_GFPp + A1_Rep5_GFPm)

#########################################################################
## Set up the model:                                                   ##
#########################################################################

## Likelihood calculation for model spread of extreme underdominance
## system through a population:

logLike_UD_Mod <- function(s, UD_GFPp, UD_GFPm) {
  
  numGens <- sum(!is.na(UD_GFPm))
  Total <- UD_GFPp + UD_GFPm
  
  ## Initial genotype numbers:
  AA_Temp <- UD_GFPp[1] ; WW_Temp <- UD_GFPm[1]
  
  W1 <- AA_Temp + WW_Temp
  GFPp_Pred <- AA_Temp / W1
  GFPm_Pred <- WW_Temp / W1

  W2 <- AA_Temp*(1-s) + WW_Temp
  AA <- AA_Temp*(1-s) / W2
  WW <- WW_Temp / W2

	for (i in 2:numGens) {

	  AA_Temp <- AA[i-1]*AA[i-1]
	  WW_Temp <- WW[i-1]*WW[i-1]
	  
	  W1 <- AA_Temp + WW_Temp
	  GFPp_Pred[i] <- AA_Temp / W1
	  GFPm_Pred[i] <- WW_Temp / W1
	  
	  W2 <- AA_Temp*(1-s) + WW_Temp
	  AA[i] <- AA_Temp*(1-s) / W2 
	  WW[i] <- WW_Temp / W2
	}

  ## Multinomial likelihood calculation:

  GFPp_Pred <- GFPp_Pred[2:numGens]
  GFPm_Pred <- GFPm_Pred[2:numGens]
  GFPp_Data <- UD_GFPp[2:numGens]
  GFPm_Data <- UD_GFPm[2:numGens]
  Total_Data <- Total[2:numGens]

  logLike <- 0
  
  for (i in 1:length(Total_Data)) {
     logLike <- (logLike + 
       GFPp_Data[i]*log(max(GFPp_Pred[i],1e-10)) + 
       GFPm_Data[i]*log(max(GFPm_Pred[i],1e-10)))
  }

  logLike
}

## Calculate log likelihood for all experiments:

logLike_AllExpts <- function(s) {

  logLikeAll <- ( logLike_UD_Mod(s, X1045A_Rep1_GFPp, X1045A_Rep1_GFPm) +
                  logLike_UD_Mod(s, X1045A_Rep2_GFPp, X1045A_Rep2_GFPm) +
                  logLike_UD_Mod(s, X1045A_Rep3_GFPp, X1045A_Rep3_GFPm) +
                  logLike_UD_Mod(s, X1045A_Rep4_GFPp, X1045A_Rep4_GFPm) +
                  logLike_UD_Mod(s, X1045A_Rep5_GFPp, X1045A_Rep5_GFPm) +
                  logLike_UD_Mod(s, X1045A_Rep6_GFPp, X1045A_Rep6_GFPm) +
                  logLike_UD_Mod(s, X1045A_Rep7_GFPp, X1045A_Rep7_GFPm) +
                  logLike_UD_Mod(s, X1045A_Rep8_GFPp, X1045A_Rep8_GFPm) +
                  logLike_UD_Mod(s, X1045A_Rep9_GFPp, X1045A_Rep9_GFPm) + 
                  logLike_UD_Mod(s, A1_Rep1_GFPp, A1_Rep1_GFPm) +
                  logLike_UD_Mod(s, A1_Rep2_GFPp, A1_Rep2_GFPm) +
                  logLike_UD_Mod(s, A1_Rep3_GFPp, A1_Rep3_GFPm) +
                  logLike_UD_Mod(s, A1_Rep4_GFPp, A1_Rep4_GFPm) +
                  logLike_UD_Mod(s, A1_Rep5_GFPp, A1_Rep5_GFPm) )

  logLikeAll
}

## Prior function:
logPrior <- function(s) {
  # Prior on s:
  logPrior_s <- dunif(s, min = 0, max = 1, log = TRUE)

  return(logPrior_s)
}

## Posterior function:
logPosterior <- function(theta) {
  ## Extract parameters from theta vector:
  s <- theta[["s"]] 
  
  ## Calculate the log prior:
  logPrior <- logPrior(s)
  
  ## Calculate the log likelihood of the data:
  logLike <- logLike_AllExpts(s)
  
  ## Calculate the posterior:
  logPosterior <- logPrior + logLike
  return(logPosterior)
}

theta <- c(s=0.5)
logPosterior(theta)

#########################################################################
## MCMC analysis:                                                      ##
#########################################################################

## Metropolis-Hastings algorithm for searching parameter space:

mcmcMH <- function(logPosterior, initTheta, proposalSD, numIterations) {
  
  # Evaluate the log posterior at initTheta:
  logPosteriorThetaCurrent <- logPosterior(initTheta)
  
  # Initialise variables:
  thetaCurrent <- initTheta
  samples <- initTheta
  accepted <- 0
  
  # Run the MCMC algorithm for numIterations interations:
  for (i in 1:numIterations) {
    # Draw a new theta from a Gaussian proposal distribution and
    # assign this to a variable called thetaProposed:
    if (proposalSD[["s"]]==0) { s_Proposed = thetaCurrent[["s"]] } else {
      s_Proposed <- rtruncnorm(n = 1, a = 0, b = 1, mean = thetaCurrent[["s"]], sd = proposalSD[["s"]])
    }
    thetaProposed <- c(s = s_Proposed)

    # Evaluate the log posterior function at the proposed theta value:
    logPosteriorThetaProposed <- logPosterior(thetaProposed)
    
    # Compute the Metropolis-Hastings (log) acceptance probability:
    logAcceptance <- logPosteriorThetaProposed - logPosteriorThetaCurrent
    
    # Use a random number to determine if thetaProposed will be accepted:
    randNum <- runif(n = 1, min = 0, max = 1)
    
    ## If accepted, update the thetaCurrent vector, etc.:
    if (randNum < exp(logAcceptance)) {
      thetaCurrent <- thetaProposed
      logPosteriorThetaCurrent <- logPosteriorThetaProposed
      accepted <- accepted + 1
    }
    
    # Add the current theta to the vector of samples:
    samples <- c(samples, thetaCurrent)
    
    # Print the current state of chain and acceptance rate:
    cat("iteration:", i, "chain:", thetaCurrent, "acceptance rate:", accepted / i, "\n")
  }
  return(samples)
}

# Running the MCMC algorithm to vary the model parameters:
mcmcTrace <- mcmcMH(logPosterior = logPosterior, # posterior distribution
                    initTheta = c(s = 0.35), # intial parameter guess
                    proposalSD = c(s = 0.001), # standard deviations of 
                    # parameters for Gaussian proposal distribution
                    numIterations = 100000) # number of iterations

trace <- matrix(mcmcTrace, ncol = 1, byrow = T)

# Use the package "coda" to convert the trace into this format:
trace <- mcmc(trace)
plot(trace)
summary(trace)

## Remove the first 10000 iterations to allow for burn-in:
traceBurn <- trace[-(1:10000),]
traceBurn <- mcmc(traceBurn)
plot(traceBurn)
summary(traceBurn)

## Check for autocorrelation:
autocorr.plot(traceBurn)

## This suggests that the chain was run for long enough for autocorrelation
## to not compromise the results.

#########################################################################
## Now let's plot the results (deterministic):                         ##
#########################################################################

## Let's use the following parameter values for the plots:
s_Fitted <- 0.348

## Trajectory calculation:

traj_UD_Mod <- function(s, UD0, numGens) {
    
   ## Initial genotype numbers:
   AA_Temp <- UD0 ; WW_Temp <- (1 - UD0)
   
   W1 <- AA_Temp + WW_Temp
   GFPp_Pred <- AA_Temp / W1
   GFPm_Pred <- WW_Temp / W1
   
   W2 <- AA_Temp*(1-s) + WW_Temp
   AA <- AA_Temp*(1-s) / W2
   WW <- WW_Temp / W2
   
   for (i in 2:numGens) {
    
     AA_Temp <- AA[i-1]*AA[i-1]
     WW_Temp <- WW[i-1]*WW[i-1]
    
     W1 <- AA_Temp + WW_Temp
     GFPp_Pred[i] <- AA_Temp / W1
     GFPm_Pred[i] <- WW_Temp / W1
     
     W2 <- AA_Temp*(1-s) + WW_Temp
     AA[i] <- AA_Temp*(1-s) / W2 
     WW[i] <- WW_Temp / W2
  }
    
  ## Prepare model predictions for output & plotting:
    
  cbind(GFPp_Pred, GFPm_Pred)
}

## Prepare the observed data for plotting (first set of plots):

GFPp_Obs_1 <- X1045A_Rep1_GFPp/Total_1045A_Rep1
GFPp_Obs_2 <- X1045A_Rep2_GFPp/Total_1045A_Rep2
GFPp_Obs_3 <- X1045A_Rep3_GFPp/Total_1045A_Rep3
GFPp_Obs_4 <- X1045A_Rep4_GFPp/Total_1045A_Rep4
GFPp_Obs_5 <- X1045A_Rep5_GFPp/Total_1045A_Rep5
GFPp_Obs_6 <- X1045A_Rep6_GFPp/Total_1045A_Rep6
GFPp_Obs_7 <- X1045A_Rep7_GFPp/Total_1045A_Rep7
GFPp_Obs_8 <- X1045A_Rep8_GFPp/Total_1045A_Rep8
GFPp_Obs_9 <- X1045A_Rep9_GFPp/Total_1045A_Rep9
GFPp_Obs_10 <- A1_Rep1_GFPp/Total_A1_Rep1
GFPp_Obs_11 <- A1_Rep2_GFPp/Total_A1_Rep2
GFPp_Obs_12 <- A1_Rep3_GFPp/Total_A1_Rep3
GFPp_Obs_13 <- A1_Rep4_GFPp/Total_A1_Rep4
GFPp_Obs_14 <- A1_Rep5_GFPp/Total_A1_Rep5

## Calculate the predicted trajectories for plotting:

UD_Pred_70 <- traj_UD_Mod(s = s_Fitted, UD0 = 0.7, numGens = 12)
UD_Pred_50 <- traj_UD_Mod(s = s_Fitted, UD0 = 0.5, numGens = 12)
UD_Pred_80 <- traj_UD_Mod(s = s_Fitted, UD0 = 0.8, numGens = 12)
UD_Pred_90 <- traj_UD_Mod(s = s_Fitted, UD0 = 0.9, numGens = 12)

GFPp_Pred_70 <- UD_Pred_70[,"GFPp_Pred"]
GFPp_Pred_50 <- UD_Pred_50[,"GFPp_Pred"]
GFPp_Pred_80 <- UD_Pred_80[,"GFPp_Pred"]
GFPp_Pred_90 <- UD_Pred_90[,"GFPp_Pred"]

## Collate observed & predicted data into a data frame:

gens <- seq(1,12)

UD_Results_GFP_70 <- data.frame(gens, GFPp_Obs_1, GFPp_Obs_2, GFPp_Obs_3, GFPp_Obs_4, 
                                GFPp_Obs_5, GFPp_Obs_6, GFPp_Pred_70)
UD_Results_GFP_50 <- data.frame(gens, GFPp_Obs_7, GFPp_Obs_8, GFPp_Obs_9, GFPp_Pred_50)
UD_Results_GFP_80 <- data.frame(gens, GFPp_Obs_10, GFPp_Obs_11, GFPp_Obs_12, GFPp_Obs_13, 
                                GFPp_Pred_50)
UD_Results_GFP_90 <- data.frame(gens, GFPp_Obs_14, GFPp_Pred_50)
UD_Results_GFP_All <- data.frame(gens, GFPp_Obs_1, GFPp_Obs_2, GFPp_Obs_3, GFPp_Obs_4, 
                                 GFPp_Obs_5, GFPp_Obs_6, GFPp_Obs_7, GFPp_Obs_8, 
                                 GFPp_Obs_9, GFPp_Obs_10, GFPp_Obs_11, GFPp_Obs_12, 
                                 GFPp_Obs_13, GFPp_Obs_14, GFPp_Pred_50, GFPp_Pred_70, 
                                 GFPp_Pred_80, GFPp_Pred_90)

## Plot the results:

ggplot(UD_Results_GFP_All, aes(x = gens, y = GFPp_Obs_1, color = Expt)) +
  geom_line(aes(y = GFPp_Obs_14, col = "90%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_10, col = "80%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_11, col = "80%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_12, col = "80%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_13, col = "80%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_1, col = "70%"), size = 1.2) + 
  geom_line(aes(y = GFPp_Obs_2, col = "70%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_3, col = "70%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_4, col = "70%"), size = 1.2) + 
  geom_line(aes(y = GFPp_Obs_5, col = "70%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_6, col = "70%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_7, col = "50%"), size = 1.2) + 
  geom_line(aes(y = GFPp_Obs_8, col = "50%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_9, col = "50%"), size = 1.2) +
  geom_line(aes(y = GFPp_Pred_90, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50, col = "50%"), linetype = "dashed", size = 0.8) +
  labs(x = "Generation", y = "Proportion of individuals GFP+") +
  ggtitle("Observed and predicted SPECIES dynamics")

#########################################################################
## Now let's plot the results (stochastic):                            ##
#########################################################################

## Let's use the following parameter values for the plots:
s_Fitted <- 0.348

sampleSize <- 50

## Trajectory calculation:

traj_UD_Mod_Stochastic <- function(s, UD0, numGens) {
  
  ## Initial genotype numbers:
  AA_Temp <- UD0 ; WW_Temp <- (1 - UD0)
  
  W1 <- AA_Temp + WW_Temp
  GFPp_Pred <- AA_Temp / W1
  GFPm_Pred <- WW_Temp / W1
  
  W2 <- AA_Temp*(1-s) + WW_Temp
  AA <- AA_Temp*(1-s) / W2
  WW <- WW_Temp / W2
  
  for (i in 2:numGens) {

    # Deterministic proportions:
    pAA <- AA[i-1]*AA[i-1] / (AA[i-1]*AA[i-1] + WW[i-1]*WW[i-1])
    pWW <- WW[i-1]*WW[i-1] / (AA[i-1]*AA[i-1] + WW[i-1]*WW[i-1])
    
    # Stochastic numbers:
    AA_Temp <- rbinom(n = 1, size = sampleSize, prob = pAA)
    WW_Temp <- sampleSize - AA_Temp
    
    # Phenotypes:
    W1 <- AA_Temp + WW_Temp
    GFPp_Pred[i] <- AA_Temp / W1
    GFPm_Pred[i] <- WW_Temp / W1
    
    # Contribution to the next generation:
    W2 <- AA_Temp*(1-s) + WW_Temp
    AA[i] <- AA_Temp*(1-s) / W2 
    WW[i] <- WW_Temp / W2
  }
  
  ## Prepare model predictions for output & plotting:
  
  cbind(GFPp_Pred, GFPm_Pred)
}

## Prepare the observed data for plotting (first set of plots):

GFPp_Obs_1 <- X1045A_Rep1_GFPp/Total_1045A_Rep1
GFPp_Obs_2 <- X1045A_Rep2_GFPp/Total_1045A_Rep2
GFPp_Obs_3 <- X1045A_Rep3_GFPp/Total_1045A_Rep3
GFPp_Obs_4 <- X1045A_Rep4_GFPp/Total_1045A_Rep4
GFPp_Obs_5 <- X1045A_Rep5_GFPp/Total_1045A_Rep5
GFPp_Obs_6 <- X1045A_Rep6_GFPp/Total_1045A_Rep6
GFPp_Obs_7 <- X1045A_Rep7_GFPp/Total_1045A_Rep7
GFPp_Obs_8 <- X1045A_Rep8_GFPp/Total_1045A_Rep8
GFPp_Obs_9 <- X1045A_Rep9_GFPp/Total_1045A_Rep9
GFPp_Obs_10 <- A1_Rep1_GFPp/Total_A1_Rep1
GFPp_Obs_11 <- A1_Rep2_GFPp/Total_A1_Rep2
GFPp_Obs_12 <- A1_Rep3_GFPp/Total_A1_Rep3
GFPp_Obs_13 <- A1_Rep4_GFPp/Total_A1_Rep4
GFPp_Obs_14 <- A1_Rep5_GFPp/Total_A1_Rep5

## Calculate the predicted trajectories for plotting:

UD_Pred_90_Rep_1 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.9, numGens = 12)
UD_Pred_90_Rep_2 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.9, numGens = 12)
UD_Pred_90_Rep_3 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.9, numGens = 12)
UD_Pred_90_Rep_4 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.9, numGens = 12)
UD_Pred_90_Rep_5 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.9, numGens = 12)
UD_Pred_90_Rep_6 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.9, numGens = 12)
UD_Pred_90_Rep_7 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.9, numGens = 12)
UD_Pred_90_Rep_8 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.9, numGens = 12)
UD_Pred_90_Rep_9 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.9, numGens = 12)
UD_Pred_90_Rep_10 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.9, numGens = 12)

UD_Pred_80_Rep_1 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.8, numGens = 12)
UD_Pred_80_Rep_2 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.8, numGens = 12)
UD_Pred_80_Rep_3 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.8, numGens = 12)
UD_Pred_80_Rep_4 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.8, numGens = 12)
UD_Pred_80_Rep_5 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.8, numGens = 12)
UD_Pred_80_Rep_6 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.8, numGens = 12)
UD_Pred_80_Rep_7 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.8, numGens = 12)
UD_Pred_80_Rep_8 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.8, numGens = 12)
UD_Pred_80_Rep_9 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.8, numGens = 12)
UD_Pred_80_Rep_10 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.8, numGens = 12)

UD_Pred_70_Rep_1 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.7, numGens = 12)
UD_Pred_70_Rep_2 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.7, numGens = 12)
UD_Pred_70_Rep_3 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.7, numGens = 12)
UD_Pred_70_Rep_4 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.7, numGens = 12)
UD_Pred_70_Rep_5 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.7, numGens = 12)
UD_Pred_70_Rep_6 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.7, numGens = 12)
UD_Pred_70_Rep_7 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.7, numGens = 12)
UD_Pred_70_Rep_8 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.7, numGens = 12)
UD_Pred_70_Rep_9 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.7, numGens = 12)
UD_Pred_70_Rep_10 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.7, numGens = 12)

UD_Pred_50_Rep_1 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.5, numGens = 12)
UD_Pred_50_Rep_2 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.5, numGens = 12)
UD_Pred_50_Rep_3 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.5, numGens = 12)
UD_Pred_50_Rep_4 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.5, numGens = 12)
UD_Pred_50_Rep_5 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.5, numGens = 12)
UD_Pred_50_Rep_6 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.5, numGens = 12)
UD_Pred_50_Rep_7 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.5, numGens = 12)
UD_Pred_50_Rep_8 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.5, numGens = 12)
UD_Pred_50_Rep_9 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.5, numGens = 12)
UD_Pred_50_Rep_10 <- traj_UD_Mod_Stochastic(s = s_Fitted, UD0 = 0.5, numGens = 12)

GFPp_Pred_90_Rep_1 <- UD_Pred_90_Rep_1[,"GFPp_Pred"]
GFPp_Pred_90_Rep_2 <- UD_Pred_90_Rep_2[,"GFPp_Pred"]
GFPp_Pred_90_Rep_3 <- UD_Pred_90_Rep_3[,"GFPp_Pred"]
GFPp_Pred_90_Rep_4 <- UD_Pred_90_Rep_4[,"GFPp_Pred"]
GFPp_Pred_90_Rep_5 <- UD_Pred_90_Rep_5[,"GFPp_Pred"]
GFPp_Pred_90_Rep_6 <- UD_Pred_90_Rep_6[,"GFPp_Pred"]
GFPp_Pred_90_Rep_7 <- UD_Pred_90_Rep_7[,"GFPp_Pred"]
GFPp_Pred_90_Rep_8 <- UD_Pred_90_Rep_8[,"GFPp_Pred"]
GFPp_Pred_90_Rep_9 <- UD_Pred_90_Rep_9[,"GFPp_Pred"]
GFPp_Pred_90_Rep_10 <- UD_Pred_90_Rep_10[,"GFPp_Pred"]

GFPp_Pred_80_Rep_1 <- UD_Pred_80_Rep_1[,"GFPp_Pred"]
GFPp_Pred_80_Rep_2 <- UD_Pred_80_Rep_2[,"GFPp_Pred"]
GFPp_Pred_80_Rep_3 <- UD_Pred_80_Rep_3[,"GFPp_Pred"]
GFPp_Pred_80_Rep_4 <- UD_Pred_80_Rep_4[,"GFPp_Pred"]
GFPp_Pred_80_Rep_5 <- UD_Pred_80_Rep_5[,"GFPp_Pred"]
GFPp_Pred_80_Rep_6 <- UD_Pred_80_Rep_6[,"GFPp_Pred"]
GFPp_Pred_80_Rep_7 <- UD_Pred_80_Rep_7[,"GFPp_Pred"]
GFPp_Pred_80_Rep_8 <- UD_Pred_80_Rep_8[,"GFPp_Pred"]
GFPp_Pred_80_Rep_9 <- UD_Pred_80_Rep_9[,"GFPp_Pred"]
GFPp_Pred_80_Rep_10 <- UD_Pred_80_Rep_10[,"GFPp_Pred"]

GFPp_Pred_70_Rep_1 <- UD_Pred_70_Rep_1[,"GFPp_Pred"]
GFPp_Pred_70_Rep_2 <- UD_Pred_70_Rep_2[,"GFPp_Pred"]
GFPp_Pred_70_Rep_3 <- UD_Pred_70_Rep_3[,"GFPp_Pred"]
GFPp_Pred_70_Rep_4 <- UD_Pred_70_Rep_4[,"GFPp_Pred"]
GFPp_Pred_70_Rep_5 <- UD_Pred_70_Rep_5[,"GFPp_Pred"]
GFPp_Pred_70_Rep_6 <- UD_Pred_70_Rep_6[,"GFPp_Pred"]
GFPp_Pred_70_Rep_7 <- UD_Pred_70_Rep_7[,"GFPp_Pred"]
GFPp_Pred_70_Rep_8 <- UD_Pred_70_Rep_8[,"GFPp_Pred"]
GFPp_Pred_70_Rep_9 <- UD_Pred_70_Rep_9[,"GFPp_Pred"]
GFPp_Pred_70_Rep_10 <- UD_Pred_70_Rep_10[,"GFPp_Pred"]

GFPp_Pred_50_Rep_1 <- UD_Pred_50_Rep_1[,"GFPp_Pred"]
GFPp_Pred_50_Rep_2 <- UD_Pred_50_Rep_2[,"GFPp_Pred"]
GFPp_Pred_50_Rep_3 <- UD_Pred_50_Rep_3[,"GFPp_Pred"]
GFPp_Pred_50_Rep_4 <- UD_Pred_50_Rep_4[,"GFPp_Pred"]
GFPp_Pred_50_Rep_5 <- UD_Pred_50_Rep_5[,"GFPp_Pred"]
GFPp_Pred_50_Rep_6 <- UD_Pred_50_Rep_6[,"GFPp_Pred"]
GFPp_Pred_50_Rep_7 <- UD_Pred_50_Rep_7[,"GFPp_Pred"]
GFPp_Pred_50_Rep_8 <- UD_Pred_50_Rep_8[,"GFPp_Pred"]
GFPp_Pred_50_Rep_9 <- UD_Pred_50_Rep_9[,"GFPp_Pred"]
GFPp_Pred_50_Rep_10 <- UD_Pred_50_Rep_10[,"GFPp_Pred"]

## Collate observed & predicted data into a data frame:

gens <- seq(1,12)

UD_Results_GFP_90 <- data.frame(gens, GFPp_Obs_14, GFPp_Pred_90_Rep_1, 
                                GFPp_Pred_90_Rep_2, GFPp_Pred_90_Rep_3, GFPp_Pred_90_Rep_4, 
                                GFPp_Pred_90_Rep_5, GFPp_Pred_90_Rep_6, GFPp_Pred_90_Rep_7, 
                                GFPp_Pred_90_Rep_8, GFPp_Pred_90_Rep_9, GFPp_Pred_90_Rep_10)
UD_Results_GFP_80 <- data.frame(gens, GFPp_Obs_10, GFPp_Obs_11, GFPp_Obs_12, GFPp_Obs_13, 
                                GFPp_Pred_80_Rep_1, 
                                GFPp_Pred_80_Rep_2, GFPp_Pred_80_Rep_3, GFPp_Pred_80_Rep_4, 
                                GFPp_Pred_80_Rep_5, GFPp_Pred_80_Rep_6, GFPp_Pred_80_Rep_7, 
                                GFPp_Pred_80_Rep_8, GFPp_Pred_80_Rep_9, GFPp_Pred_80_Rep_10)
UD_Results_GFP_70 <- data.frame(gens, GFPp_Obs_1, GFPp_Obs_2, GFPp_Obs_3, GFPp_Obs_4, 
                                GFPp_Obs_5, GFPp_Obs_6, GFPp_Pred_70_Rep_1, GFPp_Pred_70_Rep_2, 
                                GFPp_Pred_70_Rep_3, GFPp_Pred_70_Rep_4, GFPp_Pred_70_Rep_5, 
                                GFPp_Pred_70_Rep_6, GFPp_Pred_70_Rep_7, GFPp_Pred_70_Rep_8, 
                                GFPp_Pred_70_Rep_9, GFPp_Pred_70_Rep_10)
UD_Results_GFP_50 <- data.frame(gens, GFPp_Obs_7, GFPp_Obs_8, GFPp_Obs_9, GFPp_Pred_50_Rep_1, 
                                GFPp_Pred_50_Rep_2, GFPp_Pred_50_Rep_3, GFPp_Pred_50_Rep_4, 
                                GFPp_Pred_50_Rep_5, GFPp_Pred_50_Rep_6, GFPp_Pred_50_Rep_7, 
                                GFPp_Pred_50_Rep_8, GFPp_Pred_50_Rep_9, GFPp_Pred_50_Rep_10)
UD_Results_GFP_All <- data.frame(gens, GFPp_Obs_1, GFPp_Obs_2, GFPp_Obs_3, GFPp_Obs_4, 
                                 GFPp_Obs_5, GFPp_Obs_6, GFPp_Obs_7, GFPp_Obs_8, GFPp_Obs_9, 
                                 GFPp_Obs_10, GFPp_Obs_11, GFPp_Obs_12, GFPp_Obs_13, GFPp_Obs_14, 
                                 GFPp_Pred_90_Rep_1, GFPp_Pred_90_Rep_2, 
                                 GFPp_Pred_90_Rep_3, GFPp_Pred_90_Rep_4, GFPp_Pred_90_Rep_5, 
                                 GFPp_Pred_90_Rep_6, GFPp_Pred_90_Rep_7, GFPp_Pred_90_Rep_8, 
                                 GFPp_Pred_90_Rep_9, GFPp_Pred_90_Rep_10, GFPp_Pred_80_Rep_1, 
                                 GFPp_Pred_80_Rep_2, GFPp_Pred_80_Rep_3, GFPp_Pred_80_Rep_4, 
                                 GFPp_Pred_80_Rep_5, GFPp_Pred_80_Rep_6, GFPp_Pred_80_Rep_7, 
                                 GFPp_Pred_80_Rep_8, GFPp_Pred_80_Rep_9, GFPp_Pred_80_Rep_10, 
                                 GFPp_Pred_70_Rep_1, GFPp_Pred_70_Rep_2, 
                                 GFPp_Pred_70_Rep_3, GFPp_Pred_70_Rep_4, GFPp_Pred_70_Rep_5, 
                                 GFPp_Pred_70_Rep_6, GFPp_Pred_70_Rep_7, GFPp_Pred_70_Rep_8, 
                                 GFPp_Pred_70_Rep_9, GFPp_Pred_70_Rep_10, GFPp_Pred_50_Rep_1, 
                                 GFPp_Pred_50_Rep_2, GFPp_Pred_50_Rep_3, GFPp_Pred_50_Rep_4, 
                                 GFPp_Pred_50_Rep_5, GFPp_Pred_50_Rep_6, GFPp_Pred_50_Rep_7, 
                                 GFPp_Pred_50_Rep_8, GFPp_Pred_50_Rep_9, GFPp_Pred_50_Rep_10)

## Plot the results:

ggplot(UD_Results_GFP_All, aes(x = gens, y = GFPp_Obs_1, color = Expt)) +
  geom_line(aes(y = GFPp_Obs_14, col = "90%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_10, col = "80%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_11, col = "80%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_12, col = "80%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_13, col = "80%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_1, col = "70%"), size = 1.2) + 
  geom_line(aes(y = GFPp_Obs_2, col = "70%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_3, col = "70%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_4, col = "70%"), size = 1.2) + 
  geom_line(aes(y = GFPp_Obs_5, col = "70%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_6, col = "701%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_7, col = "50%"), size = 1.2) + 
  geom_line(aes(y = GFPp_Obs_8, col = "50%"), size = 1.2) +
  geom_line(aes(y = GFPp_Obs_9, col = "50%"), size = 1.2) +
  geom_line(aes(y = GFPp_Pred_90_Rep_1, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_90_Rep_2, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_90_Rep_3, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_90_Rep_4, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_90_Rep_5, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_90_Rep_6, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_90_Rep_7, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_90_Rep_8, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_90_Rep_9, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_90_Rep_10, col = "90%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80_Rep_1, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80_Rep_2, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80_Rep_3, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80_Rep_4, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80_Rep_5, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80_Rep_6, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80_Rep_7, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80_Rep_8, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80_Rep_9, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_80_Rep_10, col = "80%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70_Rep_1, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70_Rep_2, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70_Rep_3, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70_Rep_4, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70_Rep_5, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70_Rep_6, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70_Rep_7, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70_Rep_8, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70_Rep_9, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_70_Rep_10, col = "70%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50_Rep_1, col = "50%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50_Rep_2, col = "50%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50_Rep_3, col = "50%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50_Rep_4, col = "50%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50_Rep_5, col = "50%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50_Rep_6, col = "50%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50_Rep_7, col = "50%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50_Rep_8, col = "50%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50_Rep_9, col = "50%"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = GFPp_Pred_50_Rep_10, col = "50%"), linetype = "dashed", size = 0.8) +
  labs(x = "Generation", y = "Proportion of individuals GFP+") +
  ylim(0, 1) +
  scale_x_continuous(breaks=c(1, 3, 5, 7, 9, 11)) +
  ggtitle("Observed and predicted SPECIES dynamics (A1)")
