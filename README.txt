Pipeline for utilizing double digest onto a chromosome with enzymes from Bio.Restriction

Download libraries 

Make sure the libraries are insalled with:
pip install bio biopython

Create a chromosome.fasta file for double digest and choose two of your favourites enzymes. You can see all the Restriction enzymes python and a full 
documentation for dealing with them here 
http://biopython.org/DIST/docs/cookbook/Restriction.html

The main function is ddrad(chromosome, enz1, enz2). Run test.py for a first, simple example. A file e1e2.fasta will be created with all the fragments of double 
digest. Details about the length of each fragment and chromosome id are appearing in the line > bla bla bla of the fasta for meta analysis

Evolab 2023
