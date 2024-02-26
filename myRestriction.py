# cut the chromosome with a pair of enzymes and write to file.
# Everything has O(n) -linear- complexity, where n the length of chromosome

# libraries
from Bio import SeqIO, Restriction
from Bio.Seq import Seq
from Bio.Restriction import *

# Load restriction enzymes and store into RestrictionBatch
renAll = Restriction.AllEnzymes

# functions
#cut chromosome and write to file only fragments in [200,600]
# Input the chromosome.fasta file and two strings, the enzyme names
def ddrad(chromosome,e1,e2):
	enz1, enz2 = renAll.get(e1), renAll.get(e2)
	filename =  str(e1) + str(e2) + ".fasta"
	dd = ddpositions(chromosome,enz1,enz2)
	with open(filename, "w") as file1:
		for k in range(0,len(dd),2):
			if len(chromosome.seq[dd[k]-1:dd[k+1]-1]) >= 200 and len(chromosome.seq[dd[k]-1:dd[k+1]-1]) <= 600:
				file1.write(" ".join([">", str(e1), str(e2), str(k), "chromosome_id", chromosome.id, "length", str(len(chromosome.seq[dd[k]-1:dd[k+1]-1])), 
				"\n"]))
				file1.write(str(chromosome.seq[dd[k]-1:dd[k+1]-1]) + "\n")

# double digests the chromosome and ouputs the positions 
# input chromosome.fasta file and two enzymes of type Restriction in biopython 
def ddpositions(chromosome,enz1, enz2):
    enzpos1 = enz1.search(chromosome.seq)
    enzpos2 = enz2.search(chromosome.seq)
    ddpos = []; start = []; end = []
    i=0; j=0
    
    if enzpos1[0] < enzpos2[0]:
        start = enzpos1
        end = enzpos2
    else:
        start = enzpos2
        end = enzpos1

    while i < len(start) and j < len(end):
        while i < len(start) and start[i] <= end[j]:
            i+=1
        
        if i == len(start):
            ddpos.append(start[i-1])
            ddpos.append(end[j])
            return ddpos
        
        ddpos.append(start[i-1])
        ddpos.append(end[j])
        
        while  j < len(end) and  end[j] <= start[i]:
            j+=1
		
        ddpos.append(end[j-1])
        ddpos.append(start[i])
	    
        if j == len(end):
            ddpos.append(end[j-1])
            ddpos.append(start[i])
            return ddpos
