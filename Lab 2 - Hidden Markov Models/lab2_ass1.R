# install.packages("HMM")
 # library(HMM)

#Assignment 1
#(1) Build a hidden Markov model (HMM) for the scenario described above
set.seed(12345)

states = 1:10
symbols = c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")

tmatrix = matrix(0,10,10) #transmission matrix
ematrix = matrix(0,10,10) #emission matrix

tmatrix=0.5*diag(10)
ematrix=0.2*diag(10)
for(i in 1:10){
  tmatrix[i,i%%10+1]=0.5
  max=i+2
  min=i-2
  int=min:max%%10
  int[int=0]=10
  ematrix[i,int]=0.2
}

smoothingFunc <- function(hmm, obs, n){
  alpha <- exp(forward(hmm,obs))
  beta <- exp(backward(hmm,obs))
  smoothing = matrix(0,10,n)
  for(j in 1:n){
    smoothing[,j] <- beta[,j]*alpha[,j]/sum(beta[,j]*alpha[,j]) 
  }
  return(smoothing)
}

filteringFunc <- function(hmm, obs,n){
  alpha <- exp(forward(hmm,obs))
  filtering = matrix(0,10,n)
  for(j in 1:n){
    filtering[,j] <- alpha[,j]/sum(alpha[,j]) 
  }
  return(filtering)
}

hmm <- initHMM(States=states,Symbols=symbols,transProbs=tmatrix, emissionProbs=ematrix) #startprobs uniform dist.

simuHMM <- simHMM(hmm, 100)

simuHMM$states
simuHMM$observation


smoothing = smoothingFunc(hmm,simuHMM$observation,100)
filtering = filteringFunc(hmm,simuHMM$observation,100)

probPath <- viterbi(hmm,simuHMM$observation) 
smoothPath <- apply(smoothing, MARGIN=2,FUN=which.max)
filterPath <- apply(filtering, MARGIN = 2, FUN=which.max)

accuracyFunc <- function(data, inputdata){
    n=length(inputdata)
    nTrue <- sum(data == inputdata)
    return(nTrue/n)
  }

probPathAcc <- accuracyFunc(probPath,simuHMM$states)
smoothPathAcc <- accuracyFunc(smoothPath, simuHMM$states)
filterPathAcc <- accuracyFunc(filterPath, simuHMM$states)




AccSim <- function(hmm, nobs){
  sim <- simHMM(hmm, nobs)
  smooths <- smoothingFunc(hmm, sim$observation, nobs)
  filts <- filteringFunc(hmm, sim$observation, nobs)
  probpath <- viterbi(hmm, sim$observation)
  pathSmooth <- apply(smooths, MARGIN = 2, FUN=which.max)
  pathfilts <- apply(filts, MARGIN = 2, FUN=which.max)

  accpp <- accuracyFunc(probpath,sim$states)
  accfilts <- accuracyFunc(pathfilts,sim$states)
  accSmooth <- accuracyFunc(pathSmooth, sim$states)
  tempMatrix = matrix(c(accpp,accfilts,accSmooth),1,3)
  return(tempMatrix)
}

set.seed(1234566)
nSim <- 100
AccMatrix20<- matrix(0,nSim,3)
AccMatrix40 <- matrix(0,nSim,3)
AccMatrix80 <- matrix(0,nSim,3)
for(i in 1:nSim){
  AccMatrix40[i,] <- AccSim(hmm, 40)
  AccMatrix20[i,] <- AccSim(hmm,20)
  AccMatrix80[i,] <- AccSim(hmm,80)
}


plot(density(AccMatrix80[,3]) ,type="l", xlim=c(0.2,1))
lines(density(AccMatrix80[,2]), type="l")
lines(density(AccMatrix80[,1]),type="l")

# Since the smoothing takes all observations into account (backward and forward)
# the distribution has more information and is better at predicting the correct
# state. 

#

plot(density(AccMatrix80[,2]), type="l", xlim=c(0,1),)
lines(density(AccMatrix40[,2]), type="l")
lines(density(AccMatrix20[,2]), type="l")

# As seen in the plot the accuracy does not increase with the amount of observations.
# Conclusion: False. Since accuracy does not increase with amount of observations,
# amount of obs does not help us determine the next observation.

# install.packages("entropy")
# library(entropy)

#To compute step 101 simply multiply the prob. dist. for step 100 with the transition matrix
smooth101 <- t(smoothing[,100]%*%tmatrix)
smoothing[,100]

