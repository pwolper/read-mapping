import time
import os
import re                                                       #use re for easier iterate search
import csv
import sys
from datetime import datetime
from Bio import SeqIO                                           #load SeqIO from Biopython package
import matplotlib.pyplot as plt

def getGenome(file):                                    #return sequences and description from fasta file (input .fasta)                    #returns array of array
    sequences = []                                          #create empty array for for-loop
    descriptions = []                                       #create empty array for for-loop
    for record in SeqIO.parse(file, "fasta"):               #loop through all entries in source file and store sequences and descriptions
        seq = record.seq
        des = record.description
        sequences.append(str(seq))
        print(des)
        descriptions.append(des)
    return(sequences, descriptions)


def getReads(file):                                     #return reads to map from source file as array (input preferably txt)
    nucl = ["A","T","C","G"]
    reads = []
    with open(file) as read:
        for line in read.readlines():
            if any(base in line for base in nucl):
                line = line[:-1]
                reads.append(line)
    return (reads)

def compStrand(sequences):                                    #returns the complementary string of entered sequence as string
    comp_sequences = []
    nucl = {"A": "T", "C": "G", "T": "A", "G": "C"}
    for s in sequences:
        comp_s = str()
        for position,base in enumerate(s):
            if base not in nucl:
                print("Not a valid DNA sequence! Check position",position,"!")
                sys.exit("Script terminated...")
            else:
                comp_s += nucl[base]
        comp_sequences.append(comp_s[::-1])
    return comp_sequences                                     #output reverse complementary sequence



def Barplot(gotReads, mappedGenome, gotGenome, path = "./"):             #create a bar plot showing occurrence for all reads for each genome; input return from getReads, mapGenome and getGenome
    if not path =="./":                                              #create output path if not existing
        if not os.path.exists(path):
            os.makedirs(path)
    fileName = str(path)+ "/"+str(datetime.today().strftime('%Y-%m-%d_'))+ str(gotGenome[1][0]) + ".png"

    number = []
    for entry in mappedGenome:
            m = mappedGenome[entry][0] + mappedGenome[entry][1]
            number.append(len(m))
    plt.bar(range(len(number)), number)
    plt.xlabel("Read Number in File")
    plt.ylabel("Occurrences in Genome")
    plt.savefig(fileName)
    plt.clf()
    print("Boxplot created in",path,"!")

def MapSummary(mappedGenome, gotGenome, path = "./"):          #create a summary after the genome is mapped in a txt file for each genome; input return from getReads, mapGenome and getGenome
    if not path =="./":                                              #create output path if not existing
        if not os.path.exists(path):
            os.makedirs(path)
    fileName = str(path)+ "/" +str(datetime.today().strftime('%Y-%m-%d_'))+ str(gotGenome[1][0]) +".csv"            #create file name for current genome from fasta description

    with open(fileName, "w", newline = "") as f:                                     #create file
        f = csv.writer(f, delimiter = ";")
        header = ["Read", "Sense", "Antisense", "Number of sense hits", "Number of antisense hits", "Number of total hits"]
        f.writerow(header)
        for read in mappedGenome:
            sense = ",".join(str(x) for x in mappedGenome[read][0])
            anti = ",".join(str(x) for x in mappedGenome[read][1])
            n_sense = len(mappedGenome[read][0])                     #write number of found positions
            n_anti = len(mappedGenome[read][1])
            n_total = n_sense + n_anti
            cols = [str(read),sense,anti,str(n_sense),str(n_anti),str(n_total)]
            f.writerow(cols)
            # entry = ";".join(str(x) for x in cols)

    print("Successfully wrote summary file to",path,"!")

def mapGenome(genomes, Reads):                                      #return dictionary with genome description and read as key and starting positions in genome as value (returns list of list
    start = time.time()
    print("mapping reads to genome...")

    mapped = {}                                                     #create dictionary for output
    genome = str(genomes[0])
    comp_genome = str(compStrand(genomes[0]))

    for read in Reads:
          posFor = [i.start() for i in re.finditer(read, genome)]
          posRev = [i.start() for i in re.finditer(read, comp_genome)]
          posOut = [posFor, posRev]
          mapped[read] = posOut

    end = time.time()
    print("Finished mapping",len(Reads), "reads in",round(end-start),"s!")
    return(mapped)
    
    
