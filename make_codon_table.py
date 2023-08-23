# ////////////////////////////////

# This script creates a codon usage table from an input FASTA file of all
# annotated CDSs in an organism (i.e., a *_cds_from_genomic.fna file found
# at the NCBI RefSeq FTP database).

# Usage -
#       python make_codon_table.py fasta_file codon_table_filename

# Inputs -
#       1 - fasta_file = The FASTA file of CDSs from an organism
#       2 - codon_table_filename = Prefix name for the created codon table file

# Outputs - 
#       1 - A comma-delimited file of the relative frequencies (0-1)
#           of each codon from the table in the following format -
#           [Codon] [Translated amino acid] [Total number]  [Relative frequency]
#               [Codon] = triplet DNA sequence for the codon
#               [Translated amino acid] = single letter code for the amino acid
#               [Total number] =  total number of occurences for that codon
#               [Relative frequency] = codon number / total number of codons

# /////////////////////////////////

import sys

fasta_file = sys.argv[1]
table_file = sys.argv[2] + '.csv'

# Define the standard codon translation for DNA
CODONS = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
          'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
          'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
          'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
          'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
          'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
          'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
          'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
          'TAT':'Y', 'TAC':'Y', 'TAA':'X', 'TAG':'X',
          'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
          'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
          'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
          'TGT':'C', 'TGC':'C', 'TGA':'X', 'TGG':'W',
          'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
          'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
          'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

# Create trackers for coding sequences/codons not included in the codon usage analysis
# First is for partial codons, second for nonstandard stop codons, third for ambiguous nt
trackers = [0,0,0]  

def parse_codons(seq, dict, track):
    """
    Check that the input sequence is an appropriate coding sequence. 
    Then, parse the sequence for its codons and add the occurrences to a dictionary.
    """
    if len(seq)%3 != 0:
        track[0] += 1
    elif (seq[-3:] not in ['TAG', 'TAA', 'TGA']):
        track[1] += 1
    else:
        n = int(len(seq)/3)
        for i in range(n):
            codon = seq[3*i:3*i + 3]
            try:
                dict[codon] += 1
            except KeyError:
                track[2] += 1
                continue

# ////////////////////////////////////

# Create a dictionary to track counts for each codon
codon_dict = {codon:0 for codon in CODONS.keys()}
i = 0   # Tracker for the total number of DNA sequences in the input FASTA file

# Parse the input fasta file for the CDSs and count the codons in each sequence
with open(fasta_file, 'r') as ifile:
    flag = False    # Tracker to distinguish DNA sequences from title and header lines
    seq = ''
    for line in ifile:
        if line[0] == '>':
            flag = True
            i += 1
            # Check for a sequence and check that it is an appropriate CDS
            if seq != '':
                parse_codons(seq, codon_dict, trackers)
            seq = ''    # Reset the sequence for this new title line
        elif flag:
            seq = seq + line.strip().upper()
    parse_codons(seq, codon_dict, trackers)
total_number = sum(codon_dict.values())

# Write the results to an output CSV file
with open(table_file, 'w') as ofile:
    ofile.write('Codon,Amino Acid,Occurrences,Relative Frequency\n')
    for codon, num in codon_dict.items():
        ofile.write(f'{codon:s},{CODONS[codon]:s},{num:d},{num/total_number:.5f}\n')
        
# Print the total number of CDSs and number of CDSs/codons not included in the analysis
print(f'{i:d} coding sequences were parsed in the input FASTA file.')
print(f'{trackers[0]:d} coding sequences contained partial codons (not divisible by 3).')
print(f'{trackers[1]:d} coding sequences did not end with a standard stop codon.')
print(f'{trackers[2]:d} codons contained an ambiguous nucleotide(s).')