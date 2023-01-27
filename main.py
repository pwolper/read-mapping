#!/usr/bin/env python3
from time import ctime
from mapping import *

print("Currently:", time.ctime())
path = "./output"


# Running the read-mapper defined in mapping.py
genome = getGenome("data/Ecoli_genome.fasta")
Reads = getReads("data/reads_new.fasta")


# for read in Reads:
#     print([i for i in re.finditer(read,str(compStrand(genome[0])))])

mapped = mapGenome(genome,Reads)                             #mapping reads

MapSummary(mapped, genome, path)
Barplot(Reads,mapped, genome, path)
