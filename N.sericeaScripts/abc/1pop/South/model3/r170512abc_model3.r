# South population
# model3

require(abc)
n.param <- 5
prior <- read.table(file = "prior.txt", header = TRUE)
ss <- read.table(file = "sumstat_ESTSSR.txt", header = TRUE)
n.sim <- nrow(ss)
ss.data <- read.table(file = "/home/biosys/Desktop/ABC/abc/observed/1pop_South/observed_neolitsea1pop_South.txt", header = TRUE)

## range for priors
bounds <- matrix(c(1, 10000, # NCUR
                   1, 10000, # NANC
                   0, 5, # LOG10_TIME
                   0.5, 5, # MUT_SHAPE
                   0, 1), # P_MEAN
                 nrow = n.param, ncol = 2, byrow = TRUE)

## abc
res <- abc(ss.data, prior, ss,
           tol = 0.001,
           method = "neuralnet", numnet = 50,
           transf = "logit", logit.bounds = bounds)
pdf(file = "posterior.pdf")
plot(res, prior, ask = FALSE)
dev.off()

## output posterior
save(res, file = "posterior.RData")ESTSSR


