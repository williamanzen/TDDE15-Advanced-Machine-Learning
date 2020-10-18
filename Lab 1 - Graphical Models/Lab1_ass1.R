# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("RBGL")
# BiocManager::install("Rgraphviz")
# BiocManager::install("gRain")
# 
# install.packages("bnlearn")
# library(bnlearn)
# library(RBGL)
# library(Rgraphviz)
# library(gRain)

# (1)

data("asia")
hcOne <- hc(asia) #Standard HC-algorithm - comparing from this network
hcTwo <- hc(asia, start = model2network("[B][T][A|T][S|A][L|A][D|B:L][E|T:L][X|E]")) 


all.equal(hcOne,hcTwo) 
plot(hcTwo)
plot(model2network("[B][T][A|T][S|A][L|A][D|B:L][E|T:L][X|E]"))
#"Different number of directed/undirected arcs" 
#Local optimum found while using a different start model.
#Makes sense since the algorithm through iteration will end up in another solution.

hcThree <- hc(asia, restart = 1000)
all.equal(hcOne,hcThree) 
#Different arc sets found when restart = 105, 1000, otherwise no difference.
#Conclusion, there is some randomness to the HC algorithm but quite rarely does it have an impact.
#Choosing a higher restart-value will often result in a network with higher score. 

hcFour <- hc(asia, score = "bde")
all.equal(hcOne,hcFour) 
#True, i.e only changing the scoring method did not result in a different network. 
#Conclusion the scoring methods are still resulting in the highest scoring network.

hcFive <- hc(asia, score="bde", iss=5)
all.equal(hcOne,hcFive) 
#"Different number of directed/undirected arcs" iss determines the weight of the prior,
# here iss=5 compared to iss=1 which results in a different local optimum

vstructs(hcFive) #Z-Child X,Y - Parents
vstructs(hcOne)
plot(hcFive)

cpdag(hcOne) #Summary
cpdag(hcFive)
arcs(hcOne) #Printing the arcs

#(2)


#Function to compute the conditional probability of each observation
predictNet <- function(juncTree, data, features, target){
  predArray <- matrix(nrow=nrow(data),ncol=1)
  
  for(i in 1:nrow(data)){
    obsStates <- NULL
    for(p in features){
      if(data[i,p]=="yes"){
        obsStates <- c(obsStates,"yes")
      } else{
        obsStates <- c(obsStates,"no")
      }
    }
    
    
    obsEvidence <- setEvidence(object = juncTree,
                               nodes = features,
                               states = obsStates)
    
    obsPredProb <- querygrain(object = obsEvidence,
                              nodes = target)$S
    
    predArray[i] <- if(obsPredProb["yes"]>=0.5) "yes" else "no"
    
  }
  return(predArray)
}


n=dim(asia)[1]
set.seed(12345)
id=sample(1:n, floor(n*0.8))
train=asia[id,]
test=asia[-id,]

hcBNTest <- hc(train, start = model2network("[B][T][A|T][S|A][L|A][D|B:L][E|T:L][X|E]"), score="bde", iss=5)
BNTrue <- model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
plot(hcBNTest)
plot(BNTrue)
fittedTest <- bn.fit(hcBNTest, train) #Fitting parameters in structure with traindata
fittedTrue <- bn.fit(BNTrue, train)
BNGrainTest = as.grain(fittedTest) #Graphical independence network
BNGrainTrue = as.grain(fittedTrue)
junctionTreeTest <- compile(BNGrainTest) #Get cliques with potentials (Lauritzen-Spiegelhalter algorithm
predictionTrue <- predictNet(junctionTreeTrue, test, obsVars, tarVar)

confusionTest <- table(predictionTest, test$S)
confusionTest
confusionTrue <- table(predictionTrue, test$S)
confusionTrue
#Same confusion tables

# (3)


markovObsVarsTrue <- mb(BNTrue, "S")
mbPredTrue <- predictNet(junctionTreeTrue, test, markovObsVarsTrue, tarVar)

markovConfusionTrue <- table(mbPredTrue, test$S)
markovConfusionTrue
#Same confusion table as in (2) 

# (4)

naiveBayesNet <- model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")
fittedNaive <- bn.fit(naiveBayesNet, train)
naiveGrain <- as.grain(fittedNaive)
naiveJunction <- compile(naiveGrain)
mb(naiveBayesNet, "S")

predNaive <- predictNet(naiveJunction,test,obsVars,tarVar)

confusionNaive <- table(predNaive, test$S)
confusionNaive
confusionTrue

# (5)
# Markov Blanket (MB) consists of the nodes of importance for the target nodes' dependencies in the network.
# If we do not change the MB, the outcome of the prediction on the target node will not change.
# Since the MB does not change from (2) to (3) nor the two different bayesian networks used in 
# task (2), therefore the result does not change between these assignments.
# When we perform the naive bayes network, the MB has drastically changed from containing B,L to containing
# all the other nodes. Therefore we get a different result in (4).

