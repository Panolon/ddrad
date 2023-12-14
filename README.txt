Pipeline for utilizing double digest

0. Download libraries 

pip install bio biopython

1. Make a file chromosome.fasta in order to double digest it. Make sure it has only one squence

2. Choose two of your favourites enzymes. Make sure the enzymes are in python enzymes list

3. Call ddrad(chromo, e1,e2) and let the show begin!

A file e1e2.fasta will be created with all the fragments of double digest. Details about the length of each fragment 
and chromosome id are appearing in the line > bla bla bla of the fasta for meta analysis

Take care!
