#!/usr/bin/env python

"""
Directed Protein Mutagenesis Library Design Script
--------------------------

Design a custom DNA oligo library for the directed mutagenesis
of a protein's DNA sequence based on the input type of mutant(s)
to create, mutant amino acids at each residue, constant DNA oligo 
sequences, and DNA sequences to avoid. 

This script requires Python 3.7+ to be installed.

Usage -
    python mutate_protein_v1.0.py config_file.txt

Parameters in config_file.txt -
    mutate_seq     = Protein's DNA sequence to mutate
    up_const_seq   = Constant upstream DNA seqeunce to append to 
                     mutant protein's DNA sequence to create the oligo
    down_const_seq = Constant downstream DNA seqeunce to append to 
                     mutant protein's DNA sequence to create the oligo
    avoid_seqs     = DNA sequences to avoid when creating the library
    CDS_file       = FASTA file of CDSs from host organism to estimate
                     codon usage (optional)
    CUT_table      = Codon usage table for host organism (or default).
                     Ignored if CDS file given
    num_mutations  = Number of codon mutations for the library. Creates
                     libraries for any combination up to four mutations.
    name_prefix    = Prefix to append to all mutant names (optional)
    include_wt     = Yes/no to include the wildtype DNA sequence
                     in the output library
    remove_repeats = Yes/no to remove repeated DNA sequences from oligo
                     library
    output_prefix  = Prefix for output files and directory
    
    A list of residue positions and the corresponding amino acids to
    mutate at each residue that corresponds to the input DNA sequence
    

More information, including a comprehensive description of all input 
constraints, can be found in the associated README.txt file
"""

# Throughout the script, 'depth' defines the number of mutations in 
# the protein sequence. 0 = single, 1 = double, etc. 

import sys
import subprocess
import os
from shutil import move as mv
from configparser import ConfigParser
import re

config_file = sys.argv[1]
if config_file.upper() in ['H', 'HELP', '-H', '-HELP']:
    print(__doc__)
    sys.exit(0)

# Compile regular expressions to screen for appropriate inputs
unknown_nt = re.compile(r'[^ATGC]') 
unknown_aa = re.compile(r'[^FLIMVSPTAYHQNKDECWRG, ]')

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

BASE_PAIRS = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

# Default codon usage table, defined for E. coli strain MG1655 K-12
DEFAULT_CUT = {'F': [('TTT', '0.02227'), ('TTC', '0.01654')], 
               'L': [('CTG', '0.05298'), ('TTA', '0.01385'), ('TTG', '0.01364'), 
                     ('CTC', '0.01110'), ('CTT', '0.01101'), ('CTA', '0.00388')], 
               'I': [('ATT', '0.03046'), ('ATC', '0.02520'), ('ATA', '0.00423')], 
               'M': [('ATG', '0.02784')], 
               'V': [('GTG', '0.02628'), ('GTT', '0.01828'), ('GTC', '0.01529'), 
                     ('GTA', '0.01087')], 
               'S': [('AGC', '0.01607'), ('TCG', '0.00891'), ('AGT', '0.00873'), 
                     ('TCC', '0.00860'), ('TCT', '0.00841'), ('TCA', '0.00708')],
               'P': [('CCG', '0.02333'), ('CCA', '0.00843'), ('CCT', '0.00698'), 
                     ('CCC', '0.00546')], 
               'T': [('ACC', '0.02344'), ('ACG', '0.01444'), ('ACT', '0.00888'), 
                     ('ACA', '0.00699')], 
               'A': [('GCG', '0.03380'), ('GCC', '0.02559'), ('GCA', '0.02017'), 
                     ('GCT', '0.01526')], 
               'Y': [('TAT', '0.01612'), ('TAC', '0.01224')], 
               'X': [('TAA', '0.00205'), ('TGA', '0.00093'), ('TAG', '0.00023')], 
               'H': [('CAT', '0.01290'), ('CAC', '0.00970')], 
               'Q': [('CAG', '0.02894'), ('CAA', '0.01537')], 
               'N': [('AAC', '0.02162'), ('AAT', '0.01760')], 
               'K': [('AAA', '0.03368'), ('AAG', '0.01024')], 
               'D': [('GAT', '0.03218'), ('GAC', '0.01916')], 
               'E': [('GAA', '0.03964'), ('GAG', '0.01785')], 
               'C': [('TGC', '0.00644'), ('TGT', '0.00513')], 
               'W': [('TGG', '0.01526')], 
               'R': [('CGC', '0.02207'), ('CGT', '0.02099'), ('CGG', '0.00536'),
                     ('CGA', '0.00352'), ('AGA', '0.00201'), ('AGG', '0.00110')], 
               'G': [('GGC', '0.02973'), ('GGT', '0.02475'), ('GGG', '0.01105'),
                     ('GGA', '0.00786')]}

# Define different sets of amino acids
AA = list(set(CODONS.values()))
AA.remove('X')  # Remove stop codon (don't want truncated proteins)
AA_POL = ['S', 'T', 'N', 'Q']
AA_NONPOL = ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']
AA_CHARGE = ['R', 'H', 'K', 'D', 'E']
AA_NEG = ['D', 'E']
AA_POS = ['R', 'H', 'K']
AA_OTHER = ['C', 'G', 'P']

# ///////////////////////////////////////////////

def parse_config(file):
    """
    Create a dictionary of variables and list of mutations from the config file.
    """
    config = ConfigParser()
    config.read(file)
    var_dict = {}
    mutations = []
    for opt, val in config.items('parameters'):
        var_dict[opt] = val
    for residue, muts in config.items('amino acids'):
        mutations.append((residue, muts))
    return var_dict, mutations

def check_parameters(var_dict):
    """
    Check that the variable inputs are appropriate for the script
    """
    mutate_seq = True if (not unknown_nt.search(var_dict['mutate_seq'].upper())\
                          and len(var_dict['mutate_seq'])%3 == 0) else False
    up_const_seq = True if not unknown_nt.search(var_dict["up_const_seq"].upper())\
                        else False
    down_const_seq = True if not unknown_nt.search(var_dict["down_const_seq"].upper())\
                          else False
    avoid_seqs = True if not re.search(r'[^ATGCRYSWKMBDHVN ,]',\
                                       var_dict['avoid_seqs'].upper()) else False
    cds_fasta = True
    cut_table = True
    if var_dict['num_mutations'].isnumeric():
        if re.search(r'[2-9]', var_dict['num_mutations']):
            num_mutations = True if (var_dict['num_mutations'] in ['1', '2', '3', '4']) and\
                            (len(var_dict['num_mutations'].split(',')) <= len(var_dict['mutate_seq'])/3)\
                            else False
        else:
            num_mutations = True if len(var_dict['num_mutations']) <= len(var_dict['mutate_seq'])/3\
                            else False
    else: 
        if re.search(',', var_dict['num_mutations']):
            for input in var_dict['num_mutations'].split(','):
                if not input.strip() in ['1', '2', '3', '4']:
                    num_mutations = False
                    break
                num_mutations = True if len(var_dict['num_mutations'].split(',')) <= len(var_dict['mutate_seq'])/3\
                                else False
        else:
            num_mutations = False
    name_prefix = True
    include_wt = True if var_dict['include_wt'].upper() in ['Y', 'YES', 'N', 'NO']\
                      else False
    remove_repeats = True if var_dict['remove_repeats'].upper() in ['Y', 'YES', 'N', 'NO']\
                          else False
    output_prefix = True if var_dict['output_prefix'] != '' else False
    # Check the file inputs (if necesssary) by trying to open the files
    if var_dict['cds_fasta'].upper() not in ['', 'NONE', 'N']:
        try:
            open(var_dict['cds_fasta'], 'r').close()
        except FileNotFoundError:
            print('The input file for the FASTA CDS is not valid.')
            print('Check the file name in the config file.')
            sys.exit(1)
    else:
        if var_dict['cut_table'].upper() in ['DEFAULT', 'D']:
            print('Default codon usage table chosen')
        elif var_dict['cut_table'].upper() not in ['', 'NONE']:
            try:
                open(var_dict['cut_table'], 'r').close()
            except FileNotFoundError:
                print('The input file for the codon usage table is not valid.')
                print('Check the file name in the config file.')
                sys.exit(1)
        else:
            print('The default codon usage table must be specified')
            print('or a codon usage table file must be given if a CDS FASTA is not.')
            sys.exit(1)
    for key in var_dict:
        if not eval(key):
            print(f'Input {key} is incorrect. See the README file for appropriate inputs.')
            sys.exit(1)

def check_mutations(mut_lst, mut_dna_seq):
    """
    Check that the list of amino acid mutations is in the correct format 
    and is valid for the input DNA sequence.
    """
    if len(mut_lst) != len(mut_dna_seq)/3:
        print('The number of residues listed must match the number of codons')
        print('in the input DNA sequence.')
        sys.exit(1)
    keywords_lst = ['ALL', 'POLAR', 'NONPOLAR', 'CHARGED', 'NONE', 'NEGATIVE', 'POSITIVE', 'OTHER', '-']
    for i, res in enumerate(mut_lst):
        if not res[0].isnumeric():
            print(f'All residues should be a number. Please change residue {i}.')
            sys.exit(1)
        mutations = res[1].upper().strip().split(',')
        for mutation in mutations:
            if mutation.strip() not in keywords_lst:
                if unknown_aa.search(mutation):
                    print(f'The mutation(s) in residue {res[0]:s} is inappropriate.')
                    print('It must either be a list of keywords or valid amino acids.')
                    sys.exit(1)

def format_mutations(mut_lst, mut_dna_seq):
    """
    Format the list of mutations for each residue into a list of amino acids.
    """
    format_mut_lst = []
    aa_groups_dict = {'ALL': AA, 'POLAR': AA_POL, 'NONPOLAR': AA_NONPOL,
                       'CHARGED':AA_CHARGE, 'NEGATIVE':AA_NEG, 'POSITIVE':AA_POS,
                       'NONE':[], 'OTHER':AA_OTHER}
    for i, res in enumerate(mut_lst):
        muts = [mut.strip() for mut in res[1].upper().split(',')]
        temp_muts = []
        for mut in muts:
            if mut in aa_groups_dict.keys():
                temp_muts += aa_groups_dict[mut]
            elif mut == '-':
                og_aa = translate_dna(mut_dna_seq[3*i:3*i+3])
                temp_muts = AA.copy()
                temp_muts.remove(og_aa)
            else:
                temp_muts += mut
        temp_muts = list(set(temp_muts))    # Remove repeated amino acids
        format_mut_lst.append((mut_lst[i][0], temp_muts))
    return format_mut_lst

def format_avoid_seqs(avoid_seqs):
    """
    Compile a regular expression for DNA sequences to avoid 
    in the oligo library, including reverse complements.
    """
    avoid_seq_lst = [seq.strip().upper() for seq in avoid_seqs.split(',')]
    # Expand ambiguous nucleotides to all unambiguous sequences, if present
    for seq in avoid_seq_lst.copy():
        if re.search(r'[RYSWKMBDHVN]', seq):
            expanded_seqs = expand_ambiguous_nt(seq)
            avoid_seq_lst.extend(expanded_seqs)
            avoid_seq_lst.remove(seq)
    # Add the reverse complements of each sequence
    for seq in avoid_seq_lst.copy():
        avoid_seq_lst.append(get_reverse_complement(seq))
    avoid_seq_lst = list(set(avoid_seq_lst))    # Remove repeats
    re_exp = '|'.join(avoid_seq_lst)
    return re.compile(re_exp)

def expand_ambiguous_nt(dna_seq):
    """
    Expand an ambiguous nucleotide to ATGC nucleotides in a sequence.
    Return a list of all possible unambiguous DNA sequences at that 
    position or None if no ambiguous nucleotides are found.
    """
    AMBIGUOUS_NT_CODE = {'R':['A', 'G'], 'Y':['C', 'T'], 'S':['G', 'C'], 
                         'W':['A', 'T'], 'K':['G', 'T'], 'M':['A', 'C'], 
                         'B':['C', 'G', 'T'], 'D':['A', 'G', 'T'],
                         'H':['A', 'C', 'T'], 'V':['A', 'C', 'G'], 
                         'N':['A', 'G', 'C', 'T']}
    match = re.search(r'[RYSWKMBDHVN]', dna_seq)
    if match:
        num_nt = len(AMBIGUOUS_NT_CODE[match.group()])
        expanded_seqs = [dna_seq.replace(match.group(), AMBIGUOUS_NT_CODE[match.group()][i]) for i in range(num_nt)]
        # Search the expanded seqeunces for more ambiguous nucleotides
        for expand_seq in expanded_seqs.copy(): 
            if re.search(r'[RYSWKMBDHVN]', expand_seq):
                more_seqs = expand_ambiguous_nt(expand_seq)
                expanded_seqs.remove(expand_seq)
                expanded_seqs.extend(more_seqs)
        return expanded_seqs
    else:
        print('Input sequence does not have ambiguous nucleotides.')
        return None

def get_reverse_complement(dna_seq):
    """
    Return the reverse complement of the input DNA sequence.
    """
    rev_seq = ''
    for nt in dna_seq[::-1].upper():
        if nt in BASE_PAIRS.values():
            rev_seq += BASE_PAIRS[nt]
        else:
            print('Non-ATGC nucleotide detected')
            return None
    return rev_seq

def format_depth_lst(num_mutations):
    """
    Create a list of 0/1 to signify the type of mutant library to make.
    """
    if num_mutations.isnumeric() and not re.match(r'[2-9]', num_mutations):
        depth_lst = [int(depth) for depth in num_mutations]
    else:
        inputs = [input.strip().upper() for input in num_mutations.split(',')]
        depth_lst = []
        depth_lst.append(1) if ('1' in inputs)\
                            else depth_lst.append(0)
        depth_lst.append(1) if ('2' in inputs)\
                            else depth_lst.append(0)
        depth_lst.append(1) if ('3' in inputs)\
                            else depth_lst.append(0)
        depth_lst.append(1) if ('4' in inputs)\
                            else None
    return depth_lst

def import_cut_table(cut_file):
    """
    Create a dictionary of a codon usage table with the codons 
    sorted in descending order of usage.
    
    Format: {AA1:[(codon1, frequency1), (codon2, frequency2), ...], ...}
    """
    cut_dict = {}
    with open(cut_file, 'r') as cutfile:
        next(cutfile)
        for line in cutfile:
            codon = line.split(',')[0]
            aa = line.split(',')[1]
            freq = line.strip().split(',')[-1]
            if aa not in cut_dict:
                cut_dict[aa] = [(codon, freq)]
            else:
                cut_dict[aa].append((codon, freq))
    # Sort the codons for each amino acid in descending order of usage
    for aa in cut_dict:
        cut_dict[aa].sort(key = lambda x: x[1], reverse = True)
    return cut_dict
        
def count_mutants(num_muts_lst, pos1 = 0, depth = 0):
    """
    Calculate the number of mutants for the library at a given position and depth
    based on the type of library to be created.
    """
    num = 0
    # Check that mutants can be created at the given position and depth
    if pos1 <= len(num_muts_lst) - depth - 1:
        # Iterate over all remaining viable positions
        for pos in range(pos1, len(num_muts_lst) - depth):  
            if depth == 0:  # Single mutants remaining
                num += num_muts_lst[pos]
            else:   # Multiple mutants remaining
                new_depth = depth - 1
                new_pos = pos + 1
                num += num_muts_lst[pos]*count_mutants(num_muts_lst, new_pos, new_depth)
    return num

def translate_dna(dna_seq):
    """
    Translate a DNA sequence to an amino acid sequence.
    """
    # For speed, the codon dictionary is defined outside of the function 
    # and the sequence is not screened for non-ATGC characters
    if len(dna_seq)%3 == 0:
        aa_seq = ''
        n_codons = int(len(dna_seq)/3)
        for i in range(n_codons):
            codon = dna_seq[3*i:3*i+3].upper()
            aa = CODONS[codon]
            aa_seq += aa
        return aa_seq
    else:
        print('DNA sequence is not translatable.')
        return None

def rev_translate_aa_seq(aa, pos, og_dna_seq, cut_table, i):
    """
    Translate an amino acid to a DNA codon sequence using the ith
    most commonly-used codon according to a codon usage table.
    
    Return the entire translated DNA sequence. 
    
    cut_table is a dictionary in format - {AA:[codon1, codon2, ..]},
        where DNA codons are listed in descending order of frequency.
    """
    # Check that there is a codon at the ith position in the CUT table
    if i <= len(cut_table[aa]): 
        mut_codon_seq = cut_table[aa][i][0]
        dna_seq = og_dna_seq[:3*pos] + mut_codon_seq + og_dna_seq[3*(pos+1):]
        return dna_seq
    else:
        return None

def get_mutant_dna_seq(aa, pos, og_dna_seq, re_exp, cut_table, const_up_seq, const_down_seq):
    """
    Return a DNA sequence for the mutated amino acid sequence 
    based on the design criteria -
    - Try to use the most common codon for the mutated amino acid
    - Avoid the designated sequence(s) (change codon as needed)
    - If no non-conflicting DNA sequence is found when mutating 
      the desired codon, change the DNA sequence for the previous codon 
      and repeat the process
    
    Inputs - 
    aa = amino acid to reverse translate
    pos = position of AA in DNA sequence to reverse translate
    og_dna_seq = un-mutated (original) DNA sequence
    re_exp = compiled regular expression of sequences to avoid
    cut_table = codon usage table in dictionary format
    const_up_seq = constant upstream DNA sequence appended to mutated 
                   DNA sequence for screening
    const_down_seq = constant downstream DNA sequence appended to 
                     mutated DNA sequence for screening
    """
    for i in range(len(cut_table[aa])):
        # Reverse translate the DNA sequence at the corresponding position
        mut_dna_seq = rev_translate_aa_seq(aa, pos, og_dna_seq, cut_table, i)
        # Screen for sequences to avoid in input regular expression
        screen_seq = const_up_seq + mut_dna_seq + const_down_seq
        if re_exp.search(screen_seq):
            # If a sequence to avoid is found,
            if i < len(cut_table[aa])-1:
                # Check next most common codon for this amino acid
                continue
            else:
                # If no more codon sequences available at this position,
                if pos > 0:
                    # Try changing the sequence for the previous codon
                    # while using the most common codon for current position
                    new_pos = pos-1
                    new_dna_seq = og_dna_seq[:3*pos] + cut_table[aa][0][0] + og_dna_seq[3*(pos+1):]
                    new_aa = translate_dna(og_dna_seq[3*(pos-1):3*pos])
                    mut_dna_seq = get_mutant_dna_seq(new_aa, new_pos, new_dna_seq, re_exp, cut_table, const_up_seq, const_down_seq)
                    return mut_dna_seq
                else:
                    # If no possible sequence is found, return None
                    return None
        else:
            # If no sequence to avoid is found, return the sequence
            return mut_dna_seq

# Overview of this complex recursive function
# Start with mutating the first amino acid to a given amino acid
#   for that residue, avoiding the specified DNA sequences
# Check if the depth (type of mutant) indicates that the sequence 
#   should be recorded via the depth list. Length of depth list 
#   indicates maximum number of mutations for amino acid sequence
# Check the depth of mutations. If the maximum depth hasn't been reached 
#   and enough residues are remaining, recursively run the function 
#   again at a greater depth (more mutations) and position while 
#   using the newly-created mutant sequence as the input sequence
# Repeat this recursion until the full depth of the library is reached
#   at that position, then iterate over all remaining residues
# Update the dictionary of mutants if more mutant sequences were created 
#   and should be recorded
# Repeat the function for the remaining amino acids in the first residue
# Repeat the function for the remaining residues in sequential order

def make_mutants_recursive(dna_seq, mut_lst, cut_table, re_exp, const_seqs, depth = 0, depth_lst = [1,1], pos = 0, prefix = ''):
    """
    Mutate a residue in an amino acid sequence to all specified amino 
    acids and, if specified, further mutate downstream residues to 
    their respective amino acids.
    Return a dictionary of all mutant sequences in the format - 
    {Name:DNA sequence}
    
    Inputs - 
    dna_seq = input DNA sequence to mutate
    mut_lst = list of tuples containing a list of AA to mutate
        at each residue in the order of the amino acid sequence.
        Format is - [(residue1, AA_lst1), (residue2, AA_lst2), ...]
    cut_table = dictionary of codon usage in host organism in format - 
        {AA1:[(codon1, frequency1), (codon2, frequency2), ...], ...}. 
        Codons are listed in descending order of frequency
    re_exp = compiled regular expression of DNA sequence(s) to avoid
    const_seqs = list of upstream and downstream DNA sequences appended
        to mutated DNA sequence for screening of DNA sequences to avoid
    depth = numerical indicator of type of mutant being created for the
        DNA sequence
        (0 for single, 1 for double, etc.)
    depth_lst = list of 0/1 to indicate whether a mutant sequence
        should be recorded. Position corresponds to type of mutant 
        in acending order ([0,1] = record only double mutants)
    pos = amino acid residue position to start creating the mutations
    prefix = prefix to name for all mutants at the given depth
    """
    mut_dict = {}
    aa = translate_dna(dna_seq[pos*3:pos*3+3])
    num_aa = len(dna_seq)/3
    for acid in mut_lst[pos][1]:
        name = f'{prefix:s}{aa:s}{mut_lst[pos][0]:s}{acid:s}'
        mut_seq = get_mutant_dna_seq(acid, pos, dna_seq, re_exp, cut_table, *const_seqs)
        if mut_seq != None:
            # Check if type of mutant should be recorded
            if depth_lst[depth] == 1:
                mut_dict[name] = mut_seq
            # Check if higher order mutants should and can be created
            if depth != len(depth_lst)-1 and depth <= num_aa-pos:
                new_depth = depth + 1   # Increase depth
                # Iterate over all positions where creating higher mutations is possible
                for new_pos in range(pos + 1, int(num_aa)): 
                    more_muts_dict = make_mutants_recursive(mut_seq, mut_lst, cut_table, re_exp, const_seqs, new_depth, depth_lst, pos = new_pos, prefix = f'{name}-')
                    mut_dict.update(more_muts_dict)
        else:
            print(f'No mutant DNA sequence could be found for the {acid:s} amino acid at residue {mut_lst[pos][0]:d}.')
    return mut_dict 

def remove_repeats(dict):
    inv_dict = {}
    for key, val in dict.items():
        if val not in inv_dict:
           inv_dict[val] = key
    new_dict = {val:key for key, val in inv_dict.items()}
    return new_dict

def main():
    """
    Create a directed mutagenesis library for a protein region
    using the input parameters in the config file.
    """
    # Read the config file and check the inputs
    var_dict, mutations_raw = parse_config(config_file)
    check_parameters(var_dict)
    check_mutations(mutations_raw, var_dict['mutate_seq'])
    
    # Format the inputs
    mutations = format_mutations(mutations_raw, var_dict['mutate_seq'])
    re_avoid = format_avoid_seqs(var_dict['avoid_seqs'])
    depth_lst = format_depth_lst(var_dict['num_mutations'])
    mutate_dna_seq = var_dict['mutate_seq'].upper()
    
    # Define portions of the constant DNA sequences needed to screen
    # for DNA sequences to avoid, truncated in case the constant 
    # sequences contain a sequence to avoid
    max_avoid_seq_len = len(max([x.strip() for x in var_dict['avoid_seqs'].split(',')], key = len))
    const_seqs = [var_dict['up_const_seq'].upper()[-max_avoid_seq_len+1:], 
                  var_dict['down_const_seq'].upper()[:max_avoid_seq_len-1]]
    
    # Get the CUT table
    if var_dict['cds_fasta'].upper() not in ['', 'NONE']:
        # Create a codon usage table if a CDS FASTA is given
        cut_file_prefix = ''.join(var_dict['cds_fasta'].split('.')[:-1])
        print('Creating a codon usage table from the input FASTA file of CDSs...')
        subprocess.run(f'python make_codon_table.py {var_dict["cds_fasta"]} {cut_file_prefix}_cut_table')
        cut_dict = import_cut_table(f'{cut_file_prefix}_cut_table.csv')
    elif var_dict['cut_table'].upper() in ['DEFAULT', 'D']:
        print('Using the default codon usage table (E. coli MG1655)')
        cut_dict = DEFAULT_CUT
    elif var_dict['cut_table'].upper() not in ['', 'NONE']:
        print('Importing the codon usage table from the input file')
        cut_dict = import_cut_table(var_dict["cut_table"])
    
    # Estimate the total maximum size of the library 
    # (assuming all mutations can be made)       
    num_mutations = [len(mutations[i][1]) for i in range(len(mutations))]   # number mutations at each residue
    lib_size_lst = []
    for depth in range(len(depth_lst)):
        if depth_lst[depth] == 1:   # Check that mutations at this depth should be recorded
            num_muts = count_mutants(num_mutations, 0, depth)
            lib_size_lst.append(num_muts)
    if var_dict['include_wt'].upper() in ['Y', 'YES']:
        print(f'\nThe maximum library size is estimated to have {sum(lib_size_lst)+1:d} mutants.\n')
    else:
        print(f'\nThe maximum library size is estimated to have {sum(lib_size_lst):d} mutants.\n')
    
    
    # Make the library of mutant sequences by iterating over all positions
    mut_dict = {}
    if var_dict['include_wt'].upper() in ['Y', 'YES']:
        mut_dict['WT'] = var_dict['mutate_seq'].upper()
    for pos in range(len(mutations)):
        pos_mutants = make_mutants_recursive(mutate_dna_seq, mutations, cut_dict, re_avoid, const_seqs, depth = 0, depth_lst = depth_lst, pos = pos, prefix = var_dict['name_prefix'])
        mut_dict.update(pos_mutants)
    
    # Remove repeated DNA sequences from the library, if indicated
    if var_dict['remove_repeats'].upper() in ['Y', 'YES']:
        mut_dict = remove_repeats(mut_dict.copy())
        print(f'The created library has {len(mut_dict)} mutants without DNA sequence repeats.')
    else:
        print(f'The created library has {len(mut_dict)} mutants.')

    # Write the resulting mutant library to output files
    # One with details about each mutation, other the DNA oligo library
    mutation_file = var_dict['output_prefix'] + '_mutations.csv'
    library_file = var_dict['output_prefix'] + '_oligo_library.csv'
    with open(mutation_file, 'w') as mut_file,\
         open(library_file, 'w') as lib_file:
        mut_file.write('Name,DNA Sequence,AA Sequence,Original AA,Mutate AA\n')
        lib_file.write('Name,Oligo Sequence\n')
        for name in mut_dict:
            mutant = name.lstrip(var_dict['name_prefix'])
            aa_seq = translate_dna(mut_dict[name])
            mut_file.write(f'{name:s},{mut_dict[name]:s},{aa_seq:s}\n')
            oligo_seq = var_dict['up_const_seq'] + mut_dict[name] + var_dict['down_const_seq']
            lib_file.write(f'{name:s},{oligo_seq:s}\n')
    
    # Create a directory and move the files to the directory
    # If the directory already exists, keep it
    # If the files already exist, replace them
    dir_path = f'./{var_dict["output_prefix"]}'
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    for file in [mutation_file, library_file]:
        if os.path.isfile(dir_path + '/' + file):
            os.remove(dir_path + '/' + file)
        mv(file, dir_path)
    if var_dict['cds_fasta'].upper() not in ['', 'NONE']:
        cut_filename = f'{cut_file_prefix}_cut_table.csv'
        if os.path.isfile(dir_path + '/' + cut_filename):
            os.remove(dir_path + '/' + cut_filename)
        mv(cut_filename, dir_path)
    
    
# ////////////////////////////////////////////

if __name__ == '__main__':
    main()