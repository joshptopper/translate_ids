#Author: Joshua Topper

The goal of this script is to obtain the ensembl ID of the MC1R gene through
multiple queries using the restful api. The ensembl ID is the used to get the 
nucleotide sequence for the MC1R gene and write it to a fasta file. The sequence 
is used to determine the largest open reading frame (ORF) and then convert it from
DNA to amino acid. The amino acid ORF is then written to the same fasta file as the 
gene sequence. The ensembl ID will then be used to get a list of unique homologous 
species to a text file.

*NOTE: There are no input files BUT you can edit line 125 in the scriptAPI.py
script to whatever gene you require. 

To execute the script make sure you are inside the translate_ids directory and type 
"python3 scriptAPI.py" and hit enter. This will:

1) Creates a fasta file called mc1r_sequence.fasta that will contain both the original 
nucelotide sequence of the gene and the longest ORF of that gene as an amino acid sequence.

2) Creates a text file called mc1r_homology_list.txt that will contain all unique, homologous 
species of the human MC1R gene.

*NOTE: You may get this message when you execute scriptAPI.py:

/usr/local/lib64/python3.7/site-packages/Bio/Seq.py:2338: BiopythonWarning: Partial codon, len(sequence) not a multiple of three. 
Explicitly trim the sequence or add trailing N before translation. This may become an error in future.

The script will add "N"s to the end of the sequence if there is a count not divisible by three. For example, the MC1R gene
has one "N" to the end of the sequence. This was tested with and without this function and the results were the same. You can ignore 
the BiopythonWarning.