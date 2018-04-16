# split1 with full migration model
# Flat prior for LOG10_TDIV2

n.total.sim <- 1*10^6 # total number of sims (> 10^2)
n.sim.per.batch <- 100 # number of sims per batch
n.batch <- n.total.sim/n.sim.per.batch # number of batches
n.locus.ESTSSR <- 12 # number of ESTSSR loci
outfile.ESTSSR <- "sumstat_ESTSSR.txt"

## My environment
fsc <- "/home/biosys/Desktop/ABC/programs/fsc_linux64/fsc25221" # path for fastsimcoal2
tpl.ESTSSR <- "model_split1_ESTSSR.tpl" # path for tpl file of ESTSSR
arlsumstat <- "/home/biosys/Desktop/ABC/programs/arlsumstat_linux/arlsumstat3522_64bit" # path for arlsumstat
arl_run.ESTSSR <- "/home/biosys/Desktop/ABC/abc/observed/3pops/arl_run.ars" # path for arl_run.ars of ESTSSR
ssdefs.ESTSSR <- "/home/biosys/Desktop/ABC/abc/observed/3pops/ssdefs.txt" # path for ssdefs.txt of ESTSSR

## make priors (when in log scale, 5 digits are enough)
# 1) parameters for mutation model
MUT_MEAN <- 0.0001 # mean mu for ESTSSR

# shape parameter of gamma distribution related to variance of mu
MUT_SHAPE <- round(runif(n.total.sim, 0.5, 5), 3) 

# P parameter of GSM
# Upper limit of P was set to 0.5 considering the results of single population
# size change model.
P_MEAN <- round(runif(n.total.sim, 0, 0.5), 3)

# prepare variance of mu for ESTSSR
shape <- MUT_SHAPE
rate <- shape/MUT_MEAN

# prepare variance of P (Excoffier et al. 2005)
a <- 0.5+199*P_MEAN
b <- a*(1-P_MEAN)/P_MEAN

for(i in 1:n.locus.ESTSSR) {
    if(nchar(i) == 1) {
        assign(paste("P", "0", i, sep = ""), rbeta(n.total.sim, a, b))
        assign(paste("MUT", "0", i, sep = ""),
               rgamma(n.total.sim, shape, rate))
    }
    if(nchar(i) == 2) {
        assign(paste("P", i, sep = ""), rbeta(n.total.sim, a, b))
        assign(paste("MUT", i, sep = ""),
               rgamma(n.total.sim, shape, rate))
    }
}

# 2) parameters for demographic model
# divergence time (TDIV1 < TDIV2)
lower.log10.time <- 0 # lower limit of LOG10_Ts in log10(gens)
upper.log10.time <- 6 # upper limit of LOG10_Ts in log10(gens)

LOG10_TDIV2 <- round(runif(n.total.sim, lower.log10.time, upper.log10.time), 5)
R_TDIV1 <- runif(n.total.sim, 0, 1) # relative value for TDIV1 (from 0 to 1)
LOG10_TDIV1 <- round(R_TDIV1*LOG10_TDIV2, 5)

TDIV1 <- 10^LOG10_TDIV1
TDIV2 <- 10^LOG10_TDIV2

# migration rate
# Direction of migration is to coalescence (i.e. backward in time).
# MIG is a proportion of migrants per generation.
lower.log10.mig <- -5 # lower limit of LOG10_MIGs
upper.log10.mig <- -1 # lower limit of LOG10_MIGs

LOG10_MIG12 <- round(runif(n.total.sim, lower.log10.mig, upper.log10.mig), 5)
LOG10_MIG21 <- round(runif(n.total.sim, lower.log10.mig, upper.log10.mig), 5)
MIG12 <- 10^LOG10_MIG12
MIG21 <- 10^LOG10_MIG21

LOG10_MIG13 <- round(runif(n.total.sim, lower.log10.mig, upper.log10.mig), 5)
LOG10_MIG31 <- round(runif(n.total.sim, lower.log10.mig, upper.log10.mig), 5)
MIG13 <- 10^LOG10_MIG13
MIG31 <- 10^LOG10_MIG31

LOG10_MIG23 <- round(runif(n.total.sim, lower.log10.mig, upper.log10.mig), 5)
LOG10_MIG32 <- round(runif(n.total.sim, lower.log10.mig, upper.log10.mig), 5)
MIG23 <- 10^LOG10_MIG23
MIG32 <- 10^LOG10_MIG32


# 3) output priors for parameter estimation after simulation 
prior.out <- data.frame(LOG10_TDIV1, LOG10_TDIV2,
                        LOG10_MIG12, LOG10_MIG21,
                        LOG10_MIG13, LOG10_MIG31,
                        LOG10_MIG23, LOG10_MIG32,
                        MUT_SHAPE, P_MEAN)
write.table(prior.out, "prior.txt", quote = FALSE, sep = " ", row.names = FALSE)
rm(prior.out)

# 4) output priors for fsc
prior.ESTSSR <- data.frame(TDIV1, TDIV2,
                         MIG12, MIG21,
                         MIG13, MIG31,
                         MIG23, MIG32,
                         MUT01, MUT02, MUT03, MUT04, MUT05, MUT06,
                         MUT07, MUT08, MUT09, MUT10, MUT11, MUT12,
                         P01, P02, P03, P04, P05, P06,
                         P07, P08, P09, P10, P11, P12)



## simulation of genotypes and calculation for sumstas
# copy files for simulation of ESTSSR
system(paste("cp", arl_run.ESTSSR, "."))
system(paste("cp", ssdefs.ESTSSR, "."))

# simulate only one time for making header for outfile of ESTSSR
sim.dir.ESTSSR <- unlist(strsplit(tpl.ESTSSR, ".tpl"))
arp.file.ESTSSR <- paste("./", sim.dir.ESTSSR, "/", sim.dir.ESTSSR, "_1_1.arp", sep = "")
write.table(prior.ESTSSR[1, ], "temp_prior.txt",
            quote = FALSE, sep = " ", row.names = FALSE)
system(paste(fsc, "-t", tpl.ESTSSR, "-n 1 -f temp_prior.txt", "-g"))
system(paste(arlsumstat, arp.file.ESTSSR, outfile.ESTSSR, 0, 2))

# simulate per batch
# make arpfile name list
arp.file.ESTSSR <- paste("./", sim.dir.ESTSSR, "/", sim.dir.ESTSSR, "_",
                       1:n.sim.per.batch, "_1.arp", sep = "")

for(i in 1:n.batch) {
    ## get location of data in ith batch
    data.location <- ((i-1)*n.sim.per.batch+1):(i*n.sim.per.batch) 

    ## simulate ESTSSR data and calculate sumstats
    # copy files for simulation of ESTSSR
    system(paste("cp", arl_run.ESTSSR, "."))
    system(paste("cp", ssdefs.ESTSSR, "."))
    
    write.table(prior.ESTSSR[data.location, ], "temp_prior.txt",
                quote = FALSE, sep = " ", row.names = FALSE)
    system(paste(fsc, "-t", tpl.ESTSSR, "-n 1 -f temp_prior.txt", "-g"))
    for(j in 1:n.sim.per.batch) {
        system(paste(arlsumstat, arp.file.ESTSSR[j], outfile.ESTSSR, 1, 0))
    }
    # cleaning
    system(paste("rm -rf", sim.dir.ESTSSR))
    system("rm temp_prior.txt")
    system("rm *seed.txt")
    system("rm arl_run.ars")
    system("rm ssdefs.txt")
}

system("rm -f arl_pro.txt")
print("Simulations have been completed!")

