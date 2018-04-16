# model compaison among population divergence models

require(abc)

## loading observed sumstas
ss.data <- read.table(file = "/home/biosys/Desktop/ABC/abc/observed/3pops/observed_neolitsea3pops.txt", header = TRUE)

## loading simulated sumstats
ss.split1 <- read.table(file = "/home/biosys/Desktop/ABC/abc/3pops/split1/sumstat_ESTSSR.txt", header = TRUE)

ss.split1_no_mig <- read.table(file = "/home/biosys/Desktop/ABC/abc/3pops/split1_no_mig/sumstat_ESTSSR.txt", header = TRUE)

ss.split1_past_mig <- read.table(file = "/home/biosys/Desktop/ABC/abc/3pops/split1_past_mig/sumstat_ESTSSR.txt", header = TRUE)

ss.split1_mig12 <- read.table(file = "/home/biosys/Desktop/ABC/abc/3pops/split1_mig12/sumstat_ESTSSR.txt", header = TRUE)

ss.split1_past_12_mig <- read.table(file = "/home/biosys/Desktop/ABC/abc/3pops/split1_past_12_mig/sumstat_ESTSSR.txt", header = TRUE)

## preparing abc
nsim <- nrow(ss.split1)
model.index <- rep(c("split1", "split1_no_mig","split1_past_mig","split1_mig12","split1_past_12_mig"), each = nsim)
ss.model <- rbind(ss.split1, ss.split1_no_mig, ss.split1_past_mig, ss.split1_mig12,ss.split1_past_12_mig)

## abc
model.comp <- postpr(ss.data, model.index, ss.model,
                     #tol = 0.0005,
                     tol = 0.001,
                     #method = "mnlogistic")
                     method = "neuralnet", numnet = 50)
summary(model.comp)

#################################################################
# tol = 0.003
#Proportion of accepted simulations (rejection):
#            split1       split1_mig12      split1_no_mig split1_past_12_mig 
#            0.0121             0.0461             0.2059             0.1626 
#   split1_past_mig 
#            0.5733 

#Bayes factors:
#                    split1 split1_mig12 split1_no_mig split1_past_12_mig
#split1              1.0000       0.2634        0.0589             0.0746
#split1_mig12        3.7967       1.0000        0.2237             0.2833
#split1_no_mig      16.9725       4.4703        1.0000             1.2665
#split1_past_12_mig 13.4011       3.5297        0.7896             1.0000
#split1_past_mig    47.2473      12.4443        2.7837             3.5256
#                   split1_past_mig
#split1                      0.0212
#split1_mig12                0.0804
#split1_no_mig               0.3592
#split1_past_12_mig          0.2836
#split1_past_mig             1.0000


#Posterior model probabilities (neuralnet):
#            split1       split1_mig12      split1_no_mig split1_past_12_mig 
#            0.0085             0.0670             0.2373             0.2543 
#   split1_past_mig 
#            0.4328 

#Bayes factors:
#                    split1 split1_mig12 split1_no_mig split1_past_12_mig
#split1              1.0000       0.1270        0.0359             0.0335
#split1_mig12        7.8756       1.0000        0.2825             0.2636
#split1_no_mig      27.8789       3.5399        1.0000             0.9331
#split1_past_12_mig 29.8765       3.7935        1.0717             1.0000
#split1_past_mig    50.8469       6.4563        1.8238             1.7019
#                   split1_past_mig
#split1                      0.0197
#split1_mig12                0.1549
#split1_no_mig               0.5483
#split1_past_12_mig          0.5876
#split1_past_mig             1.0000

#################################################################
# tol = 0.001
#Proportion of accepted simulations (rejection):
#            split1       split1_mig12      split1_no_mig split1_past_12_mig 
#            0.0088             0.0504             0.2184             0.1798 
#   split1_past_mig 
#            0.5426 

#Bayes factors:
#                    split1 split1_mig12 split1_no_mig split1_past_12_mig
#split1              1.0000       0.1746        0.0403             0.0489
#split1_mig12        5.7273       1.0000        0.2308             0.2803
#split1_no_mig      24.8182       4.3333        1.0000             1.2147
#split1_past_12_mig 20.4318       3.5675        0.8233             1.0000
#split1_past_mig    61.6591      10.7659        2.4844             3.0178
#                   split1_past_mig
#split1                      0.0162
#split1_mig12                0.0929
#split1_no_mig               0.4025
#split1_past_12_mig          0.3314
#split1_past_mig             1.0000


#Posterior model probabilities (neuralnet):
#            split1       split1_mig12      split1_no_mig split1_past_12_mig 
#            0.0095             0.0742             0.2632             0.2500 
#   split1_past_mig 
#            0.4031 

#Bayes factors:
#                    split1 split1_mig12 split1_no_mig split1_past_12_mig
#split1              1.0000       0.1285        0.0362             0.0382
#split1_mig12        7.7834       1.0000        0.2821             0.2971
#split1_no_mig      27.5909       3.5449        1.0000             1.0530
#split1_past_12_mig 26.2017       3.3664        0.9496             1.0000
#split1_past_mig    42.2505       5.4283        1.5313             1.6125
#                   split1_past_mig
#split1                      0.0237
#split1_mig12                0.1842
#split1_no_mig               0.6530
#split1_past_12_mig          0.6201
#split1_past_mig             1.0000
#################################################################