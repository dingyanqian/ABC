# N1 population
# model comparison

require(abc)

## load observed sumstats
ss.data <- read.table(file = "/home/biosys/Desktop/ABC/abc/observed/1pop_N1/observed_neolitsea1pop_N1.txt", header = TRUE)

## load simulated sumstats
ss.model1 <- read.table(file = "/home/biosys/Desktop/ABC/abc/1pop/N1/model1/sumstat_ESTSSR.txt", header = TRUE)

ss.model2 <- read.table(file = "/home/biosys/Desktop/ABC/abc/1pop/N1/model2/sumstat_ESTSSR.txt", header = TRUE)

ss.model3 <- read.table(file = "/home/biosys/Desktop/ABC/abc/1pop/N1/model3/sumstat_ESTSSR.txt", header = TRUE)

ss.model4 <- read.table(file = "/home/biosys/Desktop/ABC/abc/1pop/N1/model4/sumstat_ESTSSR.txt", header = TRUE)

## preperation
nsim <- nrow(ss.model1)
model.index <- rep(c("model1", "model2", "model3", "model4"), each = nsim)
ss.model <- rbind(ss.model1, ss.model2, ss.model3, ss.model4)


## abc
model.comp <- postpr(ss.data, model.index, ss.model,
                     tol = 0.001,
                     method = "neuralnet", numnet = 50)
summary(model.comp)



#################################################
#tol = 0.001
#Proportion of accepted simulations (rejection):
#model1 model2 model3 model4 
#0.3990 0.0092 0.4102 0.1815 
#
#Bayes factors:
#       model1 model2 model3 model4
#model1 1.0000 43.1351  0.9726  2.1983
#model2  0.0232  1.0000  0.0225  0.0510
#model3  1.0282 44.3514  1.0000  2.2603
#model4  0.4549 19.6216  0.4424  1.0000
#
#Posterior model probabilities (neuralnet):
#model1 model2 model3 model4 
#0.4250 0.0013 0.4031 0.1705
#
#Bayes factors:
#        model1  model2  model3  model4
#model1  1.0000 317.6961   1.0545   2.4925
#model2   0.0031   1.0000   0.0033   0.0078
#model3   0.9483 301.2832   1.0000   2.3638
#model4   0.4012 127.4592   0.4231   1.0000

#################################################
#tol = 0.015
#Proportion of accepted simulations (rejection):
#model1 model2 model3 model4 
#0.3982 0.0098 0.4088 0.1832  
#
#Bayes factors:
#       model1 model2 model3 model4
#model1 1.0000 40.8462  0.9743  2.1733
#model2  0.0245  1.0000  0.0239  0.0532
#model3  1.0264 41.9231  1.0000  2.2306
#model4  0.4601 18.7949  0.4483  1.0000

#Posterior model probabilities (neuralnet):
#model1 model2 model3 model4 
#0.4148 0.0016 0.4053 0.1783
#
#Bayes factors:
#        model1  model2  model3  model4
#model1  1.0000 253.4005   1.0235   2.3263
#model2   0.0039   1.0000   0.0040   0.0092
#model3   0.9770 247.5753   1.0000   2.2729
#model4   0.4299 108.9271   0.4400   1.0000








