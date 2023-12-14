# test for double digest file

import myRestriction as mr

#read chromosome.fasta
chromo_file = "chromosome.fasta"
chromo = mr.SeqIO.read(chromo_file, "fasta")

#choose two enzymes to digest the chromosome
e1 = 'AatII'
e2 = 'EcoRI'

# cut the chromosome and write a file named 'e1e2.fasta'
mr.ddrad(chromo,e1,e2)
