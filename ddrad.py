# For analysis purpose we use the restriction enzyme package from Bio.Restriction
#
# load libraries

import pandas as pd
from Bio import SeqIO, Restriction
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq
import os
import os.path
import sys
import seaborn as sn
import pylab as pl
import numpy as np
import PyPDF2 as pdf

# Load information about restriction enzymes
renAll = Restriction.AllEnzymes
print(f"\nImported successfully {len(renAll.elements())} restriction enzymes")

# Read additional information about restriction enzymes from xlsx file
infofile = pd.read_csv("/home/loukas/Documents/bioinformatics/ddrad RE double digest.xlsx - RE characteristics.tsv", sep="\t")

print("\nRead enzymes from file\n")

palindromic_indexes = infofile["Palindromic enzymes"] == "Yes"

# Select palindromic enzymes
restr = infofile.loc[palindromic_indexes, "RE"]

print("These are the {} palindromic enzymes\n".format(len(restr)))
for e in restr:
	print(e, end="  ")

print()

# Keep all the palindromic enzymes in RestrictionBatch class

from Bio.Restriction import *

renz = RestrictionBatch([])
renz2 = []

for i in restr.index:
	if i not in [39,205]:
		renz.add(restr[i])
	else:
		print(f"\nEnzyme {restr[i]} was not found in python library")
	
        
print("\n")

# drop enzymes with same site
print("\nEnzymes with same site\n")
for i in range(len(renz.elements()) - 1):
	e1 = renz.get(renz.elements()[i])
	for j in range(i+1,len(renz.elements())):
		e2 = renz.get(renz.elements()[j])
		if e1.site == e2.site:
			print(f"{e1} : {e1.site}, {e2} : {e2.site}")
			renz2.append(e2)

for e in renz2:
	if e in renz:
		renz.remove(e)
            
print("\n\nThere are {} palindromic restriction enzymes".format(len(renz)))

# Categorize recognition sites into "frequent" if len ==4 otherwise "rare" 
frequent = RestrictionBatch([])
rare = RestrictionBatch([])

for e in renz:
	if e.size == 4:
		frequent.add(e)
	else:
		rare.add(e)

print("Enzymes frequent : {} \n{}".format(len(frequent),frequent))
print("\nEnzymes rare : {} \n{}".format(len(rare),rare))


# Load 19 chromosomes of Vitis Vinifera
fileChromo = "/home/loukas/Documents/bioinformatics/Vitis_vinifera.fa.gz"
seqs = list(SeqIO.parse(fileChromo, "fasta"))

print("\nTotal size:",sys.getsizeof(seqs))

# cut the chromosome with a pair of enzymes
# returns the positions of the double digest fragments

def double_digest(chromosome,enz1, enz2):
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
        
		if j == len(end):
			ddpos.append(end[j-1])
			ddpos.append(start[i])
			return ddpos
        
		ddpos.append(end[j-1])
		ddpos.append(start[i])
        
	return ddpos

# double digest and keep fragments with 200 <= length <= 600
def ddcut(minFragLen=200,maxFragLen=600):
	os.chdir("/home/loukas/Documents/bioinformatics")

	for ch_num in range(1,19):
		ch = seqs[ch_num]
		path = "chromosome" + str(ch.id)
		os.mkdir(path)
		os.chdir(path)
		for i in range(len(frequent.elements())):
			e1 = frequent.get(frequent.elements()[i])
			for j in range(len(rare.elements())):
				e2 = rare.get(rare.elements()[j])
				dd = double_digest(ch,e1,e2)

                #write to file
				filename = "chromosome" + ch.id + str(e1) + str(e2) + ".fasta"
				with open(filename, "w") as file1:
					for k in range(0,len(dd),2):
						if len(ch.seq[dd[k]-1:dd[k+1]-1]) >= minFragLen and len(ch.seq[dd[k]-1:dd[k+1]-1]) <= maxFragLen:
							file1.write(" ".join([">", str(e1), str(e2), str(k), "chromosome_id", ch.id, "length", str(len(ch.seq[dd[k]-1:dd[k+1]-1])), "\n"]))
							file1.write(str(ch.seq[dd[k]-1:dd[k+1]-1]) + "\n")

		os.chdir("/home/loukas/Documents/bioinformatics")

#ddcut()

# function that yields every line in fasta file without storing in RAM
import os
import os.path

def read_fasta(chromosome, enz1,enz2):
	os.chdir("/home/loukas/Documents/bioinformatics/chromosome"+str(chromosome.id))
    
	filename = "chromosome"+str(chromosome.id)+str(enz1)+str(enz2)+".fasta"
    
	if os.path.exists(filename)==False:
		filename = "chromosome"+str(chromosome.id)+str(enz2)+str(enz1)+".fasta"
    
	with open(filename, "r") as handle:
		for line in handle:
			yield line











