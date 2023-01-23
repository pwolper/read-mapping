#!/usr/bin/env python3

# Running the read-mapper defined in mapping.py
#
from mapping import *
import re

## Map the Genome and generate summary file and Barplot
# mapGenome("data/Ecoli_genome.fasta", "data/reads_new.fasta")

# genome = getGenome("data/Ecoli_genome.fasta")
genome = getGenome("data/Ecoli_genome.fasta")[0]
Gen = getGenome("data/Ecoli_genome.fasta")[1]
Reads = getReads("data/reads_new.fasta")[0:2]


print(Reads)
Comp = compStrand(Reads)
print(Comp)

mapped = {}
genome = str(genome)
# Gen = genome[1][0]

# for read in Reads:
#     posFor = [i.start() for i in re.finditer(read, genome)]
#     print("Position on forward strand:",posFor)
#     posRev = [i.start() for i in re.finditer(read, compStrand(genome))]
#     print("Position on rev strand:", posRev)
#     posOut = posFor + posRev
#     print(posOut)
#     mapped[Gen,read] = posOut

print(mapGenome(genome)))
