#Author: Joshua Topper

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from Bio.SeqUtils import six_frame_translations
from Bio.SeqUtils import nt_search

#define the base URL for the MyGene RESTful API
BASE_URL = "http://mygene.info/v3"
#define the base URL for the Ensembl REST API
ENSEMBL_REST_API_URL = "https://rest.ensembl.org"

#function to query gene by name and retrieve Entrez Gene ID
def query_gene_by_name(gene_name):
    url = f"{BASE_URL}/query?q={gene_name}&species=human"
    response = requests.get(url)
    if response.status_code == 200:
        gene_info = response.json()
        if 'hits' in gene_info and gene_info['hits']:
            return gene_info['hits'][0]['entrezgene']
    return None

#function to fetch Ensembl ID using Entrez Gene ID
def fetch_ensembl_id(entrez_gene_id):
    url = f"{BASE_URL}/gene/{entrez_gene_id}"
    response = requests.get(url)
    if response.status_code == 200:
        gene_info = response.json()
        if 'ensembl' in gene_info and 'gene' in gene_info['ensembl']:
            return gene_info['ensembl']['gene']
    return None

#function to fetch nucleotide sequence for a given Ensembl ID
def fetch_nucleotide_sequence(ensembl_id):
    url = f"{ENSEMBL_REST_API_URL}/sequence/id/{ensembl_id}?"
    response = requests.get(url, headers={"Content-Type": "text/plain"})
    if response.status_code == 200:
        nucleotide_sequence = response.text
        #check if sequence length is not divisible by 3
        if len(nucleotide_sequence) % 3 != 0:
            #calculate the number of missing bases to make it divisible by 3
            missing_bases = 3 - (len(nucleotide_sequence) % 3)
            #add trailing "N" nucleotides to make the sequence divisible by 3
            nucleotide_sequence += "N" * missing_bases
        return nucleotide_sequence
        #return response.text
    else:
        return None

#function to translate nucleotide sequence to amino acid sequence
def translate_sequence(nucleotide_sequence):
    dna_seq = Seq(nucleotide_sequence)
    amino_acid_sequence = dna_seq.translate()
    return str(amino_acid_sequence)



#function to identify start and stop codons
def find_codons(sequence):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}

    codon_indices = []
    for i in range(len(sequence) - 2):
        codon = sequence[i:i+3]
        if codon == start_codon:
            codon_indices.append(i)
        elif codon in stop_codons:
            codon_indices.append(i + 2)  #include the stop codon

    return codon_indices

#function to extract ORFs from the sequence
def extract_orfs(sequence):
    codon_indices = find_codons(sequence)
    orfs = []
    for i in range(0, len(codon_indices), 2):  #start codon to stop codon pairs
        start_index = codon_indices[i]
        stop_index = codon_indices[i+1] if i+1 < len(codon_indices) else None
        if stop_index:
            orf = sequence[start_index:stop_index + 3]
            orfs.append(orf)
    return orfs

#function to find the longest ORF
def find_longest_orf(sequence):
    orfs = extract_orfs(sequence)
    if orfs:
        return max(orfs, key=len)
    else:
        return None


#function to query homologous genes using the Ensembl REST API
def query_homologous_genes(ensembl_id):
    url = f"{ENSEMBL_REST_API_URL}/homology/id/{ensembl_id}?"
    response = requests.get(url, headers={"Content-Type": "application/json"})
    if response.status_code == 200:
        homology_data = response.json()
        homologous_species = extract_homologous_species(homology_data)
        return homologous_species
    else:
        print("Failed to fetch homology data.")
        return []

#function to extract homologous species from homology data
def extract_homologous_species(homology_data):
    homologous_species = set()
    if homology_data and 'data' in homology_data:
        for homology_entry in homology_data['data']:
            homologies = homology_entry.get('homologies', [])
            for homology in homologies:
                species = homology.get('target', {}).get('species', None)
                if species:
                    homologous_species.add(species)
    return list(homologous_species)



#main function to get Ensembl ID for MC1R gene, fetch nucleotide sequence,
#find the longest ORF, translate it, and write to a FASTA file
def main():
    gene_name = "MC1R" # <----- Change gene ID as needed
    
    #query gene by name and retrieve Entrez Gene ID
    entrez_gene_id = query_gene_by_name(gene_name)
    if entrez_gene_id:
        print(f"The Entrez Gene ID for gene {gene_name} is: {entrez_gene_id}")
        
        #fetch Ensembl ID using Entrez Gene ID
        ensembl_id = fetch_ensembl_id(entrez_gene_id)
        if ensembl_id:
            print(f"The Ensembl ID for gene {gene_name} is: {ensembl_id}")
            
            #fetch nucleotide sequence using Ensembl ID
            sequence = fetch_nucleotide_sequence(ensembl_id)
            if sequence:
                #find the longest ORF in the nucleotide sequence
                longest_orf = find_longest_orf(sequence)
                if longest_orf:
                    #print("Longest ORF:", longest_orf) #debugging
                    #translate the longest ORF to amino acid sequence
                    amino_acid_sequence = translate_sequence(longest_orf)
                    #print("Amino acid sequence:", amino_acid_sequence) #debugging
                    #write both nucleotide and amino acid sequences to a FASTA file
                    header_nucleotide = f">{gene_name} ({ensembl_id}) - Nucleotide sequence"
                    header_amino_acid = f">{gene_name} ({ensembl_id}) - Longest ORF (Amino Acid)"
                    output_file = f"{gene_name.lower()}_sequence.fasta"  #lowercase output file name
                    with open(output_file, "w") as fasta_file:
                        fasta_file.write(f"{header_nucleotide}\n{sequence}\n\n")
                        fasta_file.write(f"{header_amino_acid}\n{amino_acid_sequence}\n")
                    print(f"Nucleotide and longest ORF amino acid sequences have been written to {output_file}")
                    
                    #query homologous genes and write to file
                    homologous_species = query_homologous_genes(ensembl_id)
                    if homologous_species:
                        with open("mc1r_homology_list.txt", "w") as homology_file:
                            for species in homologous_species:
                                homology_file.write(species + "\n")
                        print("Homologous species list has been written to mc1r_homology_list.txt")
                    else:
                        print("No homologous species found.")
                else:
                    print("No ORF found in the nucleotide sequence.")
            else:
                print("Failed to fetch nucleotide sequence.")
        else:
            print(f"Ensembl ID for gene {gene_name} not found.")
    else:
        print(f"Entrez Gene ID for gene {gene_name} not found.")

if __name__ == "__main__":
    main()