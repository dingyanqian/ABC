# N2 population
# model1

require(abc)
n.param <- 3
prior <- read.table(file = "prior.txt", header = TRUE)
ss <- read.table(file = "sumstat_ESTSSR.txt", header = TRUE)
n.sim <- nrow(ss)
ss.data <- read.table(file = "/home/biosys/Desktop/ABC/abc/observed/1pop_N2/observed_neolitsea1pop_N2.txt", header = TRUE)

## range for priors
bounds <- matrix(c(1, 10000, # NCUR
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
save(res, file = "posterior.RData")


