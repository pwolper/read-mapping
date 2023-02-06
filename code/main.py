#!/usr/bin/env python3
from time import ctime
from mapping import *
import pathlib
import rpy2.robjects as robjects


print("Currently:", time.ctime())
path = "./output"


# Running the read-mapper defined in mapping.py
genome = getGenome((str(pathlib.Path().resolve())[:-4])+"/data/Ecoli_genome.fasta")
Reads = getReads((str(pathlib.Path().resolve())[:-4])+"data/reads_new.fasta")

mapped = mapGenome(genome,Reads)                             #mapping reads

file = MapSummary(mapped, genome, path)
Barplot(mapped, genome, path)



#Optional: create report


        #please enter image output path
robjects.r('''
out = "output/img"
''')
        #please enter genome.fasta input path
robjects.r('''
inputGen <- "data/Ecoli_genome.fasta"
''')
        #please enter MapSummary input path

robjects.globalenv['mapfile'] = file



robjects.r.source("Report.R")




