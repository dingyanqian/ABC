# N1 population
# model3
# Instantaneous size change model

n.total.sim <- 1*10^6 # total number of sims (> 10^2)
n.sim.per.batch <- 100 # number of sims per 1 batch
n.batch <- n.total.sim/n.sim.per.batch # number of batch
n.locus.ESTSSR <- 12 # number of loci for EST-SSR
outfile.ESTSSR <- "sumstat_ESTSSR.txt"

fsc <- "/home/biosys/Desktop/ABC/programs/fsc_linux64/fsc25221" # path for fastsimcoal2
tpl.ESTSSR <- "model3ESTSSR.tpl" # path for tpl
arlsumstat <- "/home/biosys/Desktop/ABC/programs/arlsumstat_linux/arlsumstat3522_64bit" # path for arlsumstat
arl_run <- "/home/biosys/Desktop/ABC/abc/observed/1pop_N1/arl_run.ars" # path for arl_run.ars
ssdefs <- "/home/biosys/Desktop/ABC/abc/observed/1pop_N1/ssdefs.txt" # path for ssdefs.txt

## making priors (when log scale, 5 digits are enough)
# current effective population size defined as diploid number of inds
NCUR <- as.integer(runif(n.total.sim, 1, 10000))
NANC <- as.integer(runif(n.total.sim, 1, 10000))

# time when population size changed
LOG10_TIME <- round(runif(n.total.sim, 0, 5), 5)

# shape parameter for gamma distribution related to variance of mu
MUT_SHAPE <- round(runif(n.total.sim, 0.5, 5), 3) 

# P for GSM
P_MEAN <- round(runif(n.total.sim, 0, 1), 3)

prior.out <- data.frame(NCUR, NANC, LOG10_TIME,
                        MUT_SHAPE, P_MEAN)
write.table(prior.out, "prior.txt", quote = FALSE, sep = " ", row.names = FALSE)
rm(prior.out)

## translating priors for fastsimcoal2 through -f option
NCUR_HAP <- 2*NCUR
RSIZE <- NANC/NCUR
TIME <- 10^LOG10_TIME

# adjusting variance of ESTSSR
MUT_MEAN <- 0.0001 # mean mu
shape <- MUT_SHAPE
rate <- shape/MUT_MEAN

# adjusting variance of P
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

prior.ESTSSR <- data.frame(NCUR_HAP, RSIZE, TIME,
                         MUT01, MUT02, MUT03, MUT04, MUT05, MUT06, MUT07, MUT08, MUT09, MUT10, MUT11,MUT12,
                         P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12)

## simulation and calculating sumstats
# copy files for simulation of ESTSSR
system(paste("cp", arl_run, "."))
system(paste("cp", ssdefs, "."))

# running arlsumstat one time for creating header
sim.dir.ESTSSR <- unlist(strsplit(tpl.ESTSSR, ".tpl"))
arp.file.ESTSSR <- paste("./", sim.dir.ESTSSR, "/", sim.dir.ESTSSR, "_1_1.arp", sep = "")
write.table(prior.ESTSSR[1, ], "temp_prior.txt",
            quote = FALSE, sep = " ", row.names = FALSE)
system(paste(fsc, "-t", tpl.ESTSSR, "-n 1 -f temp_prior.txt", "-g"))
system(paste(arlsumstat, arp.file.ESTSSR, outfile.ESTSSR, 0, 2))

# making arpfile name list
arp.file.ESTSSR <- paste("./", sim.dir.ESTSSR, "/", sim.dir.ESTSSR, "_",
                       1:n.sim.per.batch, "_1.arp", sep = "")

for(i in 1:n.batch) {
    data.location <- ((i-1)*n.sim.per.batch+1):(i*n.sim.per.batch)

    ## simulating ESTSSR and calculating sumstats
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
    
}

system("rm arl_run.ars")
system("rm ssdefs.txt")
system("rm -f arl_pro.txt")

print("Simulations have been completed!")
