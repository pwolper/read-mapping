#!/usr/bin/env python3

# Running the read-mapper defined in mapping.py
#

from mapping import *

## Map the Genome and generate summary file and Barplot
# mapGenome("data/Ecoli_genome.fasta", "data/reads_new.fasta")

Genome = getGenome("data/Ecoli_genome.fasta")

Reads = getReads("data/reads_new.fasta")
print(Reads[0:2])
Comp = compStrand(Reads[0:2])
print(Comp)
