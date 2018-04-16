# make input data for alrsumstat in each pop

d <- read.csv("../../data/repeat_number_data170412.csv",
              stringsAsFactors = FALSE)

locus.name <- grep("ENS", names(d), value = TRUE)
locus.name <- locus.name[seq(1, length(locus.name), by = 2)]
locus.name <- substring(locus.name, first = 1, last = nchar(locus.name)-1)

d2 <-  d[, c("Group", grep("ENS", names(d), value = TRUE))]
names(d2)[1] <- "Pop"
d2$Pop <- factor(d2$Pop,
                 levels = c("N1", "N2", "South"))
pop.name <- levels(d2$Pop)
n.pop <- nlevels(d2$Pop)
d2$Pop <- as.numeric(d2$Pop) # N1=1, N2=2, South=3
d2[is.na(d2)] <- "?" # missing data

arp1pop <- function (data, pop, outfile) {
    n.pop <- 1
    data <- data[data$Pop == pop, ]
    n.ind <- nrow(data)
    
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
    cat('SampleName="Sample 1"\n')
    cat(paste("SampleSize=", n.ind, "\n", sep = ""))
    cat("SampleData= {\n")
    
    for(i in 1:n.ind) { # roop for individual
        # 1st allele
        cat(paste("1_",  i, " 1 ", sep = ""))
        cat(unlist(data[i, paste(locus.name, "A", sep = "")]), sep = " ")
        cat("\n")
        # 2nd allele
        cat("      ")
        cat(unlist(data[i, paste(locus.name, "B", sep = "")]), sep = " ")
        cat("\n")
    }
    cat("}\n")

    cat("\n[[Structure]]\n")
    cat('StructureName="Observed data"\n')
    cat('NbGroups=1\nGroup={\n"Sample 1"\n}\n')
    sink()
}

for(i in 1:n.pop) {
    outfile <- paste("./1pop_", pop.name[i],
                     "/neolitsea1pop_", pop.name[i], ".arp",
                     sep = "")
    arp1pop(data = d2, pop = i, outfile = outfile)
}


