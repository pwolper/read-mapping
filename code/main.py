#!/usr/bin/env python3
from time import ctime
from mapping import *
import rpy2.robjects as robjects

print("Currently:", time.ctime())
path = "./output"


# Running the read-mapper defined in mapping.py
genome = getGenome("data/Ecoli_genome.fasta")
Reads = getReads("data/reads_new.fasta")

mapped = mapGenome(genome,Reads)                             #mapping reads

MapSummary(mapped, genome, path)
Barplot(Reads,mapped, genome, path)



#Run R-script


        #please enter image output path
robjects.r('''
output <- 'C:/Users/nikla/Desktop/readmap_fin/read-mapping/output/img'
''')
        #please enter genome.fasta input path
robjects.r('''
inputGenome <- "C:/Users/nikla/Desktop/readmap_fin/read-mapping/data/Ecoli_genome.fasta"
''')
        #please enter MapSummary input path
robjects.r('''
inputMap <- "C:/Users/nikla/Desktop/readmap_fin/read-mapping/output/2023-01-27_NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome.csv"
''')
robjects.r.source("read-vis.R")