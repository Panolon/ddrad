# ddrad
Code for Double Digestion Restriction Enzymes

Setting Up a Double Digest Pipeline Using Bio.Restriction Enzymes

1. Download Libraries

Ensure you have the necessary libraries installed:

pip install bio biopython

2. Create a Chromosome Sequence File

Generate a chromosome.fasta file for the double digest and select two preferred enzymes. A comprehensive list of available Restriction enzymes and detailed documentation for handling them can be found [here](http://biopython.org/DIST/docs/cookbook/Restriction.html).

3. Implement the Main Function

The primary function ddrad(chromosome, enz1, enz2) implements the double digest process. Execute test.py for a straightforward demonstration. Upon running, a file named e1e2.fasta will be generated, containing all the resulting fragments from the double digest. Each fragment's length and the corresponding chromosome ID are annotated in the fasta file for subsequent meta-analysis.

Evolab 2023
