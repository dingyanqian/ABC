#!/bin/bash
# calculate observed sumstats

# copy arlsumstat to the current directory
cp ~/Desktop/N.sericea/ABC/programs/arlsumstat_macosx/arlsumstatmac_64bit arlsumstat

infile="neolitsea1pop_South.arp"
outfile="observed_neolitsea1pop_South.txt"

./arlsumstat $infile $outfile 0 1

# remove useless files
rm ./arlsumstat
rm randseed.txt
#rm -rf *.res # if you want to the detail of calculation, comment out this line

