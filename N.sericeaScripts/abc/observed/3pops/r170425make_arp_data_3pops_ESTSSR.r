# make input data of each population for arlsumstat in nSSR
# N1=1, N2=2, South=3

d <- read.csv("repeat_number_data170412.csv",
              stringsAsFactors = FALSE)
locus.name <- grep("ENS", names(d), value = TRUE)
locus.name <- substring(locus.name, 1, nchar(locus.name)-1)
locus.name <- sort(unique(locus.name))
n.locus <- length(locus.name)

d2 <-  d[, c("Group", paste(rep(locus.name, each = 2),
                            rep(c("A", "B"), times = n.locus),
                            sep = ""))]
names(d2)[1] <- "Pop"
d2$Pop <- as.numeric(factor(d2$Pop)) # N1=1, N2=2, South=3
n.pop <- length(unique(d2$Pop))
d2[is.na(d2)] <- "?" # missing data

arp <- function (data, outfile) {
    n.pop <- length(unique(data$Pop)) # number of pops
    n.ind <- table(data$Pop)
    
    sink(outfile, append = FALSE)
    cat("[Profile]\n")
    cat('Title="Observed data"\n')
    cat(paste("NbSamples=", n.pop, "\n\n", sep = ""))
    cat("GenotypicData=1\n")
    cat("GameticPhase=0\n")
    cat("RecessiveData=0\n")
    cat("DataType=MICROSAT\n")
    cat("LocusSeparator=WHITESPACE\n")
    cat("MissingData='?'\n\n")
    
    cat("[Data]\n[[Samples]]\n")
    for (i in 1:n.pop) { # roop for populations
        cat(paste('SampleName="Sample ', i, '"\n', sep = ""))
        cat(paste("SampleSize= ", n.ind[i], "\n", sep = ""))
        cat("SampleData= {\n")
        temp <- data[data$Pop == i, ]
        for(j in 1:n.ind[i]) { # roop for individuals
            # 1st allele
            cat(paste(i, "_",  j, " 1 ", sep = ""))
            cat(unlist(temp[j, paste(locus.name, "A", sep = "")]), sep = " ")
            cat("\n")
            # 2nd allele
            cat("    ")
            cat(unlist(temp[j, paste(locus.name, "B", sep = "")]), sep = " ")
            cat("\n")
        }
        cat("}\n")
    }

    cat("\n[[Structure]]\n")
    cat('StructureName="Observed data"\n')
    cat("NbGroups=1\nGroup={\n")
    for (i in 1:n.pop) {
        cat(paste('"Sample ', i, '"\n', sep = ""))
    }
    cat("}")
    sink()
}

arp(data = d2, outfile = "neolitsea3pops.arp")




