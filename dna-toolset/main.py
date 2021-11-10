"""
Author  : Murthadza bin Aznam
Date    : 2021-11-01

This file contains to transcribe, translate, and generate proteins from a given DNA Sequence.
"""

import random
from os import path, getcwd

from DNAtoolkit import *
from utilities import color_text_on_print, read_FASTA

def main(file: str = None, FASTA: bool = False, pos: int = 0, acid_frequency: str = "M", **kwargs):
    """
    Pilih sama ada nak random string (Tiada input), 
        atau nak baca dari FASTA (Hanya boleh output satu sahaja), 
        atau dari fail selain fasta (macam .txt) yang hanya mengandungi string DNA.

    :params acid_frequency: pilih asid amino mana yang nak dikira
    """
    length = kwargs.get("length")
    if file == None:
        rndDNAStr = ''.join([random.choice(DNA_Bases) for nuc in range(length)])
        DNAStr = validate_sequence(rndDNAStr)
    elif FASTA == True:
        file = path.join(getcwd(), 'dataset', file)
        FASTAData = read_FASTA(file)
        DNAList = []
        for i in list(FASTAData.values()):
            DNAList.append(validate_sequence(i))
        DNAStr = DNAList[pos]
    else:
        file = path.join(getcwd(), 'dataset', file)
        with open(file, 'r') as f:
            DNAStr = f.read()

    print(f"\nSequence: {color_text_on_print(DNAStr)}\n")
    print(f"[1] + Sequence length: {len(DNAStr)}\n")

    print(color_text_on_print(f"[2] + Nucleotide Frequency: {count_base_frequency(DNAStr)}\n"))

    print(f"[3] + DNA/RNA Transcription: {color_text_on_print(make_transcription(DNAStr))}\n")

    print(f"[4] + DNA String + Reverse Complement:\n5' {color_text_on_print(DNAStr)} 3'")
    print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
    print(f"3' {color_text_on_print(make_complement(DNAStr))} 5'  [Complement]")
    print(f"5' {color_text_on_print(make_reverse_complement(DNAStr))} 3'  [Rev. Complement]\n")

    print(f"[5] + GC Content: {calculate_gc_content(DNAStr)}%\n")
    print(f"[6] + GC Content in Subsection k=5: {calculate_gc_content_of_subsections(DNAStr, k=5)}\n")

    print(f"[7] + Amino Acids Sequence from DNA: {make_translation(DNAStr, 0)}\n")

    select_acid = acid_frequency
    print(f"[8] + Codon Frequency ({select_acid}): {calculate_codon_frequency(DNAStr, select_acid)}\n")

    print("[9] + Reading frames:")
    for frame in make_reading_frames(DNAStr):
        print("\t", frame)

    print("\n[10] + All proteins in 6 open reading frames:")
    for prot in make_proteins_from_open_reading_frames(DNAStr, 0, 0, True):
        print("\t", prot)

if __name__ == "__main__":
    main(acid_frequency="_", length=50)