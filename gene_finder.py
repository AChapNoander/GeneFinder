# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna_master = load_seq("./data/X73525.fa")


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
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
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
    sequence
    dna: a DNA sequence represented as a string
    returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
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
    list_dna.append(dna)            # Adds the remaining characters back in
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
    >>> find_all_ORFs_oneframe("GCATGAATGTAG")
    ['ATG']
    """
    list_dna = []
    hold = dna
    length = int(len(dna)/3)
    for i in range(0, length):      # Seperates out string into codons
        list_dna.append(dna[:3])    # Adds codon to list
        dna = dna[3:]               # Modifies original DNA string
    list_dna.append(dna)            # Adds the remaining characters back in
    dna = hold
    indexes = []
    for i in range(0, len(list_dna)):
        if(list_dna[i] == 'ATG'):
            indexes.append(i*3)
    to_return = []
    hold_index = -1
    for i in range(0, len(indexes)):
        if(hold_index < indexes[i]):
            dna = dna[indexes[i]:]
            to_append = rest_of_ORF(dna)
            to_return.append(to_append)
            hold_index = indexes[i] + len(to_append)
        else:
            break
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
    solution = []
    first = find_all_ORFs_oneframe(dna)
    second = find_all_ORFs_oneframe(dna[1:])
    third = find_all_ORFs_oneframe(dna[2:])
    to_add = [first, second, third]
    solution = []
    for i in to_add:
        if i is not None:
            solution += i
    return solution


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    base_pair = get_reverse_complement(dna)
    base = find_all_ORFs(dna)
    paired = find_all_ORFs(base_pair)
    solution = []
    to_add = [base, paired]
    solution = []
    for i in to_add:
        if i is not None:
            solution += i
    return solution


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
    longest = 0
    to_test = []
    for i in range(0, num_trials):
        shuffled_dna = shuffle_string(dna)
        to_test += longest_ORF(shuffled_dna)
        longest = len(max(to_test))
    return longest


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
    solution = ''
    for i in range(0, int(len(dna)/3)):
        tested = dna[:3]
        dna = dna[3:]
        solution += aa_table[tested]
    return solution


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    longest = find_all_ORFs_both_strands(dna)
    print('THRESHOLD:', threshold)
    threshold = 600
    protiens = []
    for i in longest:
        print('LEN:', len(i))
        if(len(i) > threshold):
            protiens.append(coding_strand_to_AA(i))
    print(protiens)
    return protiens


if __name__ == "__main__":
    gene_finder(dna_master)
    import doctest
    doctest.testmod()
