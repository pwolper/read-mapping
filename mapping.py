import time
from datetime import datetime
from Bio import SeqIO                                           #load SeqIO from Biopython package
import matplotlib.pyplot as plt
import os
import re                                                       #use re for easier iterate search

def getGenome(file):                                    #return sequences and description from fasta file (input .fasta)                    #returns array of array
    sequences = []                                          #create empty array for for-loop
    descriptions = []                                       #create empty array for for-loop
    for record in SeqIO.parse(file, "fasta"):               #loop through all entries in source file and store sequences and descriptions
        seq = record.seq
        des = record.description
        sequences.append(str(seq))
        descriptions.append(des)
    return(sequences, descriptions)


def getReads(file):                                     #return reads to map from source file as array (input preferably txt)
    reads = []
    with open(file) as read:
        for line in read.readlines():
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
                return "Not a valid DNA sequence!"
            else:
                comp_s += nucl[base]
        comp_sequences.append(comp_s[::-1])
    return comp_sequences                                     #output reverse complementary sequence



def Barplot(gotReads, mappedGenome, gotGenome, path):             #create a bar plot showing occurrence for all reads for each genome; input return from getReads, mapGenome and getGenome
    number = []
    # l = len(gotGenome[1])
    # for i in range(l):
    for entry in gotReads:
            m = mappedGenome[entry][0] + mappedGenome[entry][1]
            number.append(len(m))
    fileName = str(path)+ "/"+str(datetime.today().strftime('%Y-%m-%d_'))+ str(gotGenome[1]) + ".png"
    plt.bar(range(len(number)), number)
    plt.xlabel("Read Number in File")
    plt.ylabel("Occurrences in Genome")
    plt.savefig(fileName)
    plt.clf()
    print("Boxplot created!")

def MapSummary(gotReads, mappedGenome, gotGenome, path):          #create a summary after the genome is mapped in a txt file for each genome; input return from getReads, mapGenome and getGenome
    summary = {}
    # l = len(gotGenome[1])
    # for i in range(l):
    fileName = str(path)+ "/" +str(datetime.today().strftime('%Y-%m-%d_'))+ str(gotGenome[1]) +".csv"            #create file name for current genome from fasta description
    for Read in gotReads:
            summary[Read] = mappedGenome[Read]      #create dictionary for read and positions (array)
            with open(fileName, "w") as f:                                     #create file
                f.write("Read; Forward Start 5' -> 3'; Reverse Start 5' -> 3'; Number of Hits \n")
                for element in summary:
                    f.write(element)                                        #write read
                    f.write("; ")
                    f.write(", ".join(str(x) for x in summary[element][0]))    #write forward positions, separated by ";"
                    f.write("; ")
                    f.write(", ".join(str(x) for x in summary[element][1]))    #write backward positions, separated by ";"
                    f.write("; [")
                    f.write(str(len(summary[element])))                     #write number of found positions
                    f.write("]")
                    f.write("\n")
    print("Successfully wrote summary file!")

def mapGenome(genomes, Reads, path ="./"):                                      #return dictionary with genome description and read as key and starting positions in genome as value (returns list of list
    start = time.time()
    if not path =="./":                                              #create output path if not existing
        if not os.path.exists(path):
            os.makedirs(path)
    # Genomes = getGenome(G_file)                                     #store sequences and description from fasta files (input .fasta)
    # print("Finished evaluating 'getGenomes'!")
    # Reads = getReads(R_file)                                        #store reads to map from source file as array (input preferably txt)
    # print("Finished evaluating 'getReads'!")
    mapped = {}                                                     #create dictionary for output
    # g = 0                                                           #
    # for genome in genomes[0]:
    Gen = genomes[1]
    #    g = +1
    genome = str(genomes[0])
      #  compgenome = compStrand(genome)
    for read in Reads:
          posFor = [i.start() for i in re.finditer(read, genome)]
          posRev = [i.start() for i in re.finditer(read, compStrand(genome))]
          posOut = [posFor, posRev]
          mapped[read] = posOut
    # Barplot(Reads,mapped, Genomes, path)
    # MapSummary(Reads,mapped, Genomes, path)
    end = time.time()
    print("Finished mapping reads in",round(end-start),"s!")
    return(mapped)
    
    
