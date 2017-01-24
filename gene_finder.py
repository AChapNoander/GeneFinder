# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return 'NA'


def get_reverse_complement(dna):
    reversed_dna = dna[::-1]            # reverses string
    list_dna = list(reversed_dna)       # creates list from string
    blank_to_return = ''
    for i in range(0, len(list_dna)):
        blank_to_return += get_complement(list_dna[i])
    return blank_to_return


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    stop_codons = ['TAG', 'TAA', 'TGA']
    list_dna = []
    hold = dna
    length = int(len(dna)/3)
    for i in range(0, length):      # Seperates out string into codons
        list_dna.append(dna[:3])    # Adds codon to list
        dna = dna[3:]               # Modifies original DNA string
    list_dna.append(dna[::])        # Adds the remaining characters back in
    final_index = -1                # Establishes base case
    for i in range(0, len(list_dna)):
        for o in range(0, 3):
            if(list_dna[i] == stop_codons[o]):
                final_index = i*3
    dna = hold                      # Done to allow for full printing w/o stop
    if(final_index == -1):          # Done to fix dropping of last character
        return dna
    else:
        return dna[:final_index]


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    to_return = []
    index = 0
    while(index > -1):                          # Until no more start codons
        index = dna.find('ATG')                 # Searches for start codon
        if(index == -1):
            break
        if(index % 3 != 0):
            dna = dna[index+(3-index % 3):]
            index = dna.find('ATG')
            if(index == -1 or index % 3 != 0):
                break
            else:
                dna = dna[index:]
        else:
            dna = dna[index:]
        to_append = rest_of_ORF(dna)
        to_return.append(to_append)
        length = len(to_append)
        dna = dna[length:]
    return to_return


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    a = 1
    solution = []
    try:
        solution.append(find_all_ORFs_oneframe(dna)[0])
    except IndexError:
        a = 2
    try:
        solution.append(find_all_ORFs_oneframe(dna[1:])[0])
    except IndexError:
        a = 3
    try:
        solution.append(find_all_ORFs_oneframe(dna[2:])[0])
    except IndexError:
        a = 4
    return solution
    print(a)


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    base_pair = get_reverse_complement(dna)
    solution = []
    a = 1
    try:
        solution.append(find_all_ORFs(dna)[0])
    except IndexError:
        a = 2
    try:
        solution.append(find_all_ORFs(base_pair)[0])
    except IndexError:
        a = 3
    return solution
    print(a)


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    solution = find_all_ORFs_both_strands(dna)
    longest = ''
    for i in range(0, len(solution)):
        if(len(solution[i]) > len(longest)):
            longest = solution[i]
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    I = ['I', 'ATT', 'ATC', 'ATA']
    L = ['L', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG']
    V = ['V', 'GTT', 'GTC', 'GTA', 'GTG']
    F = ['F', 'TTT', 'TTC']
    M = ['M', 'ATG']
    C = ['C', 'TGT', 'TGC']
    A = ['A', 'GCT', 'GCC', 'GCA', 'GCG']
    G = ['G', 'GGT', 'GGC', 'GGA', 'GGG']
    P = ['P', 'CCT', 'CCC', 'CCA', 'CCG']
    T = ['T', 'ACT', 'ACC', 'ACA', 'ACG']
    S = ['S', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']
    Y = ['TAT', 'TAC']
    W = ['W', 'TGG']
    Q = ['Q', 'CAA', 'CAG']
    N = ['N', 'AAT', 'AAC']
    H = ['H', 'CAT', 'CAC']
    E = ['E', 'GAA', 'GAG']
    D = ['D', 'GAT', 'GAC']
    K = ['K', 'AAA', 'AAG']
    R = ['R', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
    amino = [I, L, V, F, M, C, A, G, P, T, S, Y, W, Q, N, H, E, D, K, R]
    solution = ''
    for i in range(0, int(len(dna)/3)):
        tested = dna[:3]
        dna = dna[3:]
        for o in range(0, len(amino)):
            for q in range(0, len(amino[o])):
                if(tested == amino[o][q]):
                    solution += amino[o][0]
    return solution


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()
