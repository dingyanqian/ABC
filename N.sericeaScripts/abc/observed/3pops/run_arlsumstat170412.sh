#!/bin/bash
# calculate observed sumstats

# copy arlsumstat to the current directory
cp ~/Desktop/N.sericea/ABC/programs/arlsumstat_macosx/arlsumstatmac_64bit arlsumstat

infile="neolitsea3pops.arp"
outfile="observed_neolitsea3pops.txt"

./arlsumstat $infile $outfile 0 1

# remove useless files
rm ./arlsumstat
rm randseed.txt
#rm -rf *.res # if you want to the detail of calculation, comment out this line

