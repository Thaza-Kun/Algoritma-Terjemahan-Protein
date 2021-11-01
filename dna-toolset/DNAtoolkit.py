"""
Author  : Murthadza bin Aznam
Date    : 2021-11-01

This file contains functions to transcribe, translate, and generate proteins from a given DNA Sequence.
"""
from typing import Counter
from structures import *

def validateSeq(dna_seq: str) -> str:
    """
    Check the sequence to make sure it is a DNA String
    :return: str OR False
    """
    temp_seq = dna_seq.upper()
    for nuc in temp_seq:
        if nuc not in Nucleotides:
            return False
    return temp_seq

def countNucFreq(dna_seq: str) -> dict:
    """
    Counts the frequency of each nucleotides in dna string
    :return: dict
    """
    tempFreqDict = {
        "A": 0,
        "C": 0,
        "G": 0,
        "T": 0
    }
    for nuc in dna_seq:
        tempFreqDict[nuc] += 1
    return tempFreqDict

def transcription(dna_seq: str) -> str:
    """
    Replaces thymine with uracil
    """
    return dna_seq.replace('T','U')

def complement(dna_seq: str) -> str:
    """
    Generates a complement dna string
    """
    # General Approach
    # return ''.join([DNA_ReverseComplement[nuc] for nuc in dna_seq])[::-1]

    # Pythonic Approach
    mapping = str.maketrans('ATGC', 'TACG')
    return dna_seq.translate(mapping)

def reverse_complement(dna_seq: str) -> str:
    """
    Generates a reverse complement dna string
    """
    return complement(dna_seq)[::-1]

def gc_content(dna_seq: str, decimal: int = 2) -> float:
    """
    Calculates the percentages of G and C in a given DNA string
    """
    return round((dna_seq.count('C') + dna_seq.count('G')) / len(dna_seq) * 100, decimal)

def gc_content_subsec(dna_seq: str, k=20) -> list:
    """
    Calculates the percentages of G and C in chunks of DNA strings
    :return: list of gc content percentages
    """
    res = []
    for i in range(0, len(dna_seq) - k + 1, k):
        subseq = dna_seq[i: i + k]
        res.append(gc_content(subseq))
    return res

def translate_seq(dna_seq: str, init_pos: int = 0) -> list:
    """
    Translate a given DNA string into a list of amino acids from DNA_Codon mapping
    :return: list of amino acids
    """
    res = list()
    for pos in range(init_pos, len(dna_seq) - 2, 3):
        res.append(DNA_Codons[dna_seq[pos:pos + 3]])
    return res

def codon_usage(dna_seq: str, amino_acid: str) -> dict:
    """
    Calculates the frequency (in percentages) of each codon that generates a given amino acid
    :return: dict of codon frequency
    """
    templist = list()
    for i in range(0, len(dna_seq) - 2, 3):
        if DNA_Codons[dna_seq[i:i + 3]] == amino_acid:
            templist.append(dna_seq[i:i + 3])

    freqDict = dict(Counter(templist))
    totalWight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWight, 2)
    return freqDict

def gen_reading_frames(dna_seq: str):
    """
    Generates the six reading frames from a DNA sequence, including the reverse complement
    """
    frames = []
    frames.append(translate_seq(dna_seq, 0))
    frames.append(translate_seq(dna_seq, 1))
    frames.append(translate_seq(dna_seq, 2))
    frames.append(translate_seq(reverse_complement(dna_seq), 0))
    frames.append(translate_seq(reverse_complement(dna_seq), 1))
    frames.append(translate_seq(reverse_complement(dna_seq), 2))

    return frames

def proteins_from_rf(amino_acid_seq: list) -> list:
    """
    Compute all possible proteins in an amino acid seq and return a list of all possible proteins
    :return: list of proteins
    """
    current_protein = []
    proteins = []
    for amino_acid in amino_acid_seq:
        if amino_acid == "_":
            # STOP Codon
            if current_protein:
                # Only when there is a current protein
                for p in current_protein:
                    proteins.append(p)
                current_protein = []
        else:
            if amino_acid == "M":
                # START Codon
                current_protein.append("") #? This thing adds length to the list, so that the for loop can start only after M.
            for i in range(len(current_protein)):
                # Accumulate everything else
                current_protein[i] += amino_acid
    return proteins

def all_proteins_from_orfs(dna_seq, start_read_pos=0, end_read_pos=0, ordered=False):
    """
    --- Compute all possible proteins for all open reading frames
    --- Protein Search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2
    --- API can be used to pull protein info
    """
    if end_read_pos > start_read_pos:
        rfs = gen_reading_frames(dna_seq[start_read_pos: end_read_pos])
    else:
        rfs = gen_reading_frames(dna_seq)

    res = []
    for rf in rfs:
        proteins = proteins_from_rf(rf)
        for p in proteins:
            res.append(p)

    if ordered:
        res = sorted(res, key=len, reverse=True)
        
    return res