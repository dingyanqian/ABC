# split1 model

require(abc)
prior <- read.table(file = "prior.txt", header = TRUE)
n.param <- ncol(prior)
ss <- read.table(file = "sumstat_ESTSSR.txt", header = TRUE)
n.sim <- nrow(ss)
ss.data <- read.table(file = "/home/biosys/Desktop/ABC/abc/observed/3pops/observed_neolitsea3pops.txt", header = TRUE)

## range for priors
n.time <- 2
lower.log10.time <- 0
upper.log10.time <- 6
n.mig <- 4
lower.log10.mig <- -5
upper.log10.mig <- -1


bounds <- matrix(c(rep(c(lower.log10.time, upper.log10.time), # LOG10_TIMEs
                       times = n.time),
                   rep(c(lower.log10.mig, upper.log10.mig), # LOG10_MIGs
                       times = n.mig),
                   0.5, 5, # MUT_SHAPE
                   0, 0.5), # P_MEAN
                 nrow = n.param, ncol = 2, byrow = TRUE)

## abc
res <- abc(ss.data, prior, ss,
           #tol = 0.004,
           tol = 0.001,
           method = "neuralnet", numnet = 50,
           transf = "logit", logit.bounds = bounds)
pdf(file = "posterior.pdf")
plot(res, prior, ask = FALSE)
dev.off()

## output posterior for drawing figures
save(res, file = "posterior.RData")

