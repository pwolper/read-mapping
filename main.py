#!/usr/bin/env python3

# Running the read-mapper defined in mapping.py
#
from mapping import *

path = "./output"

genome = getGenome("data/Ecoli_genome.fasta")
Reads = getReads("data/reads_new.fasta")[0:50]

mapped = mapGenome(genome,Reads,path)                             #mapping reads

MapSummary(Reads, mapped, genome, path)
Barplot(Reads,mapped, genome, path)
