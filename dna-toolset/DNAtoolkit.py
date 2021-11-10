"""
Author  : Murthadza bin Aznam
Date    : 2021-11-01

This file contains functions to transcribe, translate, and generate proteins from a given DNA Sequence.
"""
from typing import Counter
from structures import *

def validate_sequence(dna_seq: str) -> str:
    """
    Check the sequence to make sure it is a DNA String
    :return: str OR False
    """
    temp_seq = dna_seq.upper()
    for nuc in temp_seq:
        if nuc not in DNA_Bases:
            return False
    return temp_seq

def count_base_frequency(dna_seq: str) -> dict:
    """
    Counts the frequency of each nucleotides in dna string
    :return: dict
    """
    # bakul = bekas yang akan simpan sejumlah objek yang sepadan dengan label
    bakul = {
        "A": 0,
        "C": 0,
        "G": 0,
        "T": 0
    }
    for nuc in dna_seq:
        bakul[nuc] += 1
    return bakul

def make_transcription(dna_seq: str) -> str:
    """
    Replaces thymine with uracil
    """
    return dna_seq.replace('T','U')

def make_complement(dna_seq: str) -> str:
    """
    Generates a complement dna string
    """
    # General Approach
    # return ''.join([DNA_ReverseComplement[nuc] for nuc in dna_seq])[::-1]

    # Pythonic Approach
    mapping = str.maketrans('ATGC', 'TACG')
    return dna_seq.translate(mapping)

def make_reverse_complement(dna_seq: str) -> str:
    """
    Generates a reverse complement dna string
    """
    return make_complement(dna_seq)[::-1]

def calculate_gc_content(dna_seq: str, decimal: int = 2) -> float:
    """
    Calculates the percentages of G and C in a given DNA string
    """
    return round((dna_seq.count('C') + dna_seq.count('G')) / len(dna_seq) * 100, decimal)

def calculate_gc_content_of_subsections(dna_seq: str, k=20) -> list:
    """
    Calculates the percentages of G and C in chunks of DNA strings
    :return: list of gc content percentages
    """
    res = []
    for i in range(0, len(dna_seq) - k + 1, k):
        subseq = dna_seq[i: i + k]
        res.append(calculate_gc_content(subseq))
    return res

def make_translation(dna_seq: str, init_pos: int = 0) -> list:
    """
    Translate a given DNA string into a list of amino acids from DNA_Codon mapping
    :return: list of amino acids
    """
    res = list()
    for pos in range(init_pos, len(dna_seq) - 2, 3):
        res.append(DNA_Codons_Mapping[dna_seq[pos:pos + 3]])
    return res

def calculate_codon_frequency(dna_seq: str, amino_acid: str) -> dict:
    """
    Calculates the frequency (in percentages) of each codon that generates a given amino acid
    :return: dict of codon frequency
    """
    templist = list()
    for i in range(0, len(dna_seq) - 2, 3):
        if DNA_Codons_Mapping[dna_seq[i:i + 3]] == amino_acid:
            templist.append(dna_seq[i:i + 3])

    freqDict = dict(Counter(templist))
    totalWight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWight, 2)
    return freqDict

def make_reading_frames(dna_seq: str):
    """
    Generates the six reading frames from a DNA sequence, including the reverse complement
    """
    frames = []
    frames.append(make_translation(dna_seq, 0))
    frames.append(make_translation(dna_seq, 1))
    frames.append(make_translation(dna_seq, 2))
    frames.append(make_translation(make_reverse_complement(dna_seq), 0))
    frames.append(make_translation(make_reverse_complement(dna_seq), 1))
    frames.append(make_translation(make_reverse_complement(dna_seq), 2))

    return frames

def make_proteins_from_reading_frame(amino_acid_seq: list) -> list:
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

def make_proteins_from_open_reading_frames(dna_seq, start_read_pos=0, end_read_pos=0, ordered=False):
    """
    --- Compute all possible proteins for all open reading frames
    --- Protein Search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2
    --- API can be used to pull protein info
    """
    if end_read_pos > start_read_pos:
        rfs = make_reading_frames(dna_seq[start_read_pos: end_read_pos])
    else:
        rfs = make_reading_frames(dna_seq)

    res = []
    for rf in rfs:
        proteins = make_proteins_from_reading_frame(rf)
        for p in proteins:
            res.append(p)

    if ordered:
        res = sorted(res, key=len, reverse=True)
        
    return res