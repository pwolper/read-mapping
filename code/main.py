#!/usr/bin/env python3
from time import ctime
from mapping import *
import pathlib
import rpy2.robjects as robjects



path = str(pathlib.Path().resolve().parents[0]) + "/output"                               #enter output path
genomePath = "/data/Ecoli_genome.fasta"         #enter reference genome location
readPath = "/data/reads_new.fasta"               #enter read location

createReport = False                             #enter boolean value, Do you want to create analysis and report?
repOut = "output/img"                           #enter report image outputfolder


print("Currently:", time.ctime())

# Running the read-mapper defined in mapping.py
genome = getGenome((str(pathlib.Path().resolve().parents[0]))+genomePath)
Reads = getReads((str(pathlib.Path().resolve().parents[0]))+readPath)

mapped = mapGenome(genome,Reads)                             #mapping reads

file = MapSummary(mapped, genome, path)
Barplot(mapped, genome, path)



#Optional: create report

if not createReport:
    pass
else:
    robjects.globalenv['out'] = repOut
    robjects.globalenv['inputGen'] = genomePath
    robjects.globalenv['mapfile'] = file
    robjects.r.source("Report.R")




