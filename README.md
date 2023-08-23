# **Directed protein mutagenesis oligo library design script**

This script designs DNA oligonucleotide pool libraries for custom directed mutagenesis libraries for a portion of a protein sequence. The script mutates the input residues of a protein DNA sequence to the specified amino acids at each residue, using the most common codon for that amino acid for the host organism, if possible. Created sequences are screened for input enzyme recognition sequence(s) that may be used for DNA assembly of the oligo pool into a destination vector (or other DNA sequences that should be avoided). The script will output the DNA and amino acid sequences for each designed protein mutant and append constant upstream and downstream DNA sequences to create the oligo pool library. An additional option to remove repeated DNA sequences is provided, if desired. This script uses the standard amino acid translation table and can create mutagenesis libraries with up to quadruple mutants. 


--------------------
## Requirements to run the script

1. Python version 3.7 or higher

Configuration files to recreate the libraries used for this work are provided and can be ran to ensure the script is working properly.

Runtime for custom libraries will vary depending on the input conditions, with higher order libraries and larger protein residue sequences taking longer. The example library above should run relatively quickly (within 5 seconds on a Windows 10 computer with an Intel i5-10300H CPU and 8 GB RAM). 


--------------------
## Running the script

The script is run by placing all input files and the scripts in the same directory and using a command line to run the script. The following command should be typed on the command line to start the script - 
```
python make_directed_mutagenesis_library_v1.0.py config_file.txt
```
where "config_file.txt" is the name of the configuration file with your custom parameters. Note that for some systems, the command to call Python may differ slightly (i.e., Mac users may need to use "python3").


--------------------
## Inputs and outputs

Parameters to create custom libraries are specified in the configuration file and include - 

- mutate_seq = DNA sequence of the protein amino acids to mutate. Must be ATCG characters, either case, and the length must be divisible by 3
- up_const_seq = Constant upstream DNA sequence for the DNA oligo to append to the mutated protein DNA sequence for the oligo pool. Must be ATGC characters
- down_const_seq = Constant downstream DNA sequence for the DNA oligo to append to the mutated protein DNA sequence for the oligo pool. Must be ATGC characters
- avoid_seqs = DNA sequences (i.e., enzyme recognition site) that should be avoided in the creation of the oligo pool library. The script automatically includes the reverse complement. Can be any - nucleotide symbols (including ambiguous), with multiple inputs separated by commas
- CDS_file = FASTA file name of coding sequences for the host organism to predict codon usage. If present, the codon usage will be predicted from the complete CDSs in the input file using the accompanying make_codon_table.py script. If blank or "none", the script will check the input for the CUT_table parameter
- CUT_table = Indication of what codon usage table to use. Either a file name for a CSV file of codon usage for the host organism (see example file for the correct format) or "default"/"d" for the default codon usage table (E. coli MG1655 K-12) coded in the script
- num_mutations = Number of codon mutations to create for the library. Any combination of 1, 2, 3, or 4, with multiple inputs separated by a comma. Also accepts inputs as a series of 1s (include) and 0s (exclude), with the first position indicating 1 mutation, second 2 mutations, etc. For example, "101" indicates to create a library of mutants with 1 and 3 mutations only
- name_prefix = Prefix for all mutant names. Can be blank or any alphanumeric character accepted by Python
- include_wt = Indication to include the wildtype protein DNA sequence in the library. Either "yes" or "no" (or their single letters)
- remove_repeats = Indication to remove repeated DNA sequences from the library. Either "yes" or "no" (or their single letters)
- output_prefix = Prefix to output file names and name of the output directory

The residue positions and their corresponding amino acids to mutate to are input below the parameters. Residue positions must be numbers, and the number of residues must match the number of codons input for mutate_seq. The amino acids for each residue are either an amino acid single letter code or keyword, with multiple inputs separated by commas. Keywords and the related amino acid are as follows:
- all = All 20 common amino acids (F, L, I, M, V, S, P, T, A, Y, H, Q, N, K, D, E, C, W, R, G)
- none = Do not mutate at this residue. Can also be blank
- polar = Amino acids with polar R groups (S, T, N, Q)
- nonpolar = Amino acids with nonpolar R groups (A, I, L, M, F, W, Y, V)
- charged = Amino acids with charged R groups (R, H, K, D, E)
- positive = Amino acids wtih positively charged (basic) R groups (R, H, K)
- negative = Amino acids with negatively charged (acidic) R groups (D, E)
- other = Amino acids with R groups that do not fit into these other categories (C, G, P)
- \- = All 19 amino acids that are not native to this residue (i.e., if G is the amino acid at the residue, all amino acids except G)


Depending on the input parameters, the script will output two or three files in a directory with the same name as the output_prefix parameter. Two files are included for all configuration files - a CSV file (\*_mutations.csv) containing details on each mutant created in the library (name, DNA and amino acid sequences) and a CSV file (\*_oligo_library.csv) containing the corresponding DNA oligo sequence for each mutant for an oligo pool to create the directed protein mutagenesis library.


--------------------
## Script procedure - 

The procedure for the main script is as follows, with a corresponding detailed flowchart given in the Supplemental Information of the manuscript - 

1. The configuration file is parsed and input parameters validated according to the above criteria. If a parameter is invalid, the script will print the invalid parameter to the console and exit. The user should review the input parameter constraints and adjust as needed. If all input parameters are valid, the mutations, type of library, and sequences to avoid are formated appropriately. 

2. The codon usage table is imported or created. If a FASTA file of CDSs is given in the CDS_file parameter, the codon usage table is created (*_cut_table.csv) using the supplementary make_codon_table.py script. If a CSV file is given in the CUT_table parameter, the data is imported. If the default table is indicated, it is used.

3. The maximum library size, based on the number of amino acid mutations at each residue (assuming no repeated DNA sequence) and number of mutation(s) that should be created, is calculated and printed to the console. This also accounts for the wildtype sequence, if it should be included.

4. The library of mutants is created using a recursive function that iterates over every residue position and corresponding list of mutant amino acids for each residue to ensure full coverage. For each mutant sequence created, the script checks if the order of mutant should be recorded and records the DNA sequence, if needed. The process works as follows - 
   1. Beginning at the first residue, mutate to the first mutant amino acid. 
   2. If higher order mutants are indicated and enough residues are left to create higher order mutants, the residue in the next position is mutated to its first mutant amino acid. 
   3. This process continues until the highest order mutant is created, then the residue at this position is mutated to all corresponding mutant amino acids. 
   4. This residue is reverted to its original sequence, and the next residue is then mutated to all its corresponding mutant amino acids. This process continues until the last residue in the amino acid sequence is mutated to its corresponding mutant amino acids. 
   5. The script will then return to the previous lower order ("base") mutant at its specified residue, and mutate this residue to its next mutant amino acid. 
   6. The same process for creating higher order mutants is repeated. 
   7. The process of iterating through each residue, mutating to one of its corresponding mutant amino acids, then mutating over the next residue(s) to create higher order mutants (if enough residues are left) is repeated until all residues are mutated as indicated.

5. If indicated, repeated DNA sequences are removed from the final library. The size of the final library of mutants is then printed to the console.

6. The constant upstream and downstream DNA sequences to create the oligo library are appended to each mutant's DNA sequence. The DNA oligo pool is written to a CSV file (\*_oligo_library.csv), and the library of mutants and their key details (name, DNA sequence, amino acid sequence) are written to another CSV file (\*_mutations.csv). An output directory is created and these files transferred here. If a codon usage table was created from a FASTA file of CDSs, this file is also transferred to the output directory.


For the supplementary script to calculated codon usage (make_codon_table.py) - 

1. Each DNA sequence in the FASTA file is inspected to determine if it is a valid protein coding sequence. For the purposes of this script, this means that the DNA sequence is composed of complete codons (i.e., the total length is divisible by 3) and the sequence ends with a valid stop codon (TAG, TAA, TGA). Any DNA sequence that does not meet these two criteria is discarded when determining codon usage.

2. If the DNA sequence meets the two criteria above, each codon in the sequence is translated to the corresponding amino acid sequence based on the standard amino acid translation table. The total number of occurances for each codon for each amino acid is recorded. If the codon contains any ambiguous nucleotides (R, Y, S, W, K, M, B, D, H, V, or N), then the codon is ignored.

3. After iterating over all DNA sequences, the details and results for each codon (amino acid translation, total number of occurances, and ratio of occurance of the codon to total number of codons) are written to an output CSV file. 

4. The total number DNA sequences parsed and number of DNA sequences and codons removed from analysis are printed to the console.
