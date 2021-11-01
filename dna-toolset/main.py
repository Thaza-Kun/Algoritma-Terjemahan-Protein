"""
Author  : Murthadza bin Aznam
Date    : 2021-11-01

This file contains to transcribe, translate, and generate proteins from a given DNA Sequence.
"""

import random
from os import path, getcwd

from DNAtoolkit import *
from utilities import colored, readFASTA

def main(file: str = None, FASTA: bool = False, pos: int = 0):
    """
    Pilih sama ada nak random string (Tiada input), 
        atau nak baca dari FASTA (Hanya boleh output satu sahaja), 
        atau dari fail selain fasta (macam .txt) yang hanya mengandungi string DNA.
    """
    if file == None:
        rndDNAStr = ''.join([random.choice(Nucleotides) for nuc in range(50)])
        DNAStr = validateSeq(rndDNAStr)
    elif FASTA == True:
        file = path.join(getcwd(), 'dataset', file)
        FASTAData = readFASTA(file)
        DNAList = []
        for i in list(FASTAData.values()):
            DNAList.append(validateSeq(i))
        DNAStr = DNAList[pos]
    else:
        file = path.join(getcwd(), 'dataset', file)
        with open(file, 'r') as f:
            DNAStr = f.read()

    print(f"\nSequence: {colored(DNAStr)}\n")
    print(f"[1] + Sequence length: {len(DNAStr)}\n")

    print(colored(f"[2] + Nucleotide Frequency: {countNucFreq(DNAStr)}\n"))

    print(f"[3] + DNA/RNA Transcription: {colored(transcription(DNAStr))}\n")

    print(f"[4] + DNA String + Reverse Complement:\n5' {colored(DNAStr)} 3'")
    print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
    print(f"3' {colored(complement(DNAStr))} 5'  [Complement]")
    print(f"5' {colored(reverse_complement(DNAStr))} 3'  [Rev. Complement]\n")

    print(f"[5] + GC Content: {gc_content(DNAStr)}%\n")
    print(f"[6] + GC Content in Subsection k=5: {gc_content_subsec(DNAStr, k=5)}\n")

    print(f"[7] + Amino Acids Sequence from DNA: {translate_seq(DNAStr, 0)}\n")

    select_acid = "M"
    print(f"[8] + Codon Frequency ({select_acid}): {codon_usage(DNAStr, select_acid)}\n")

    print("[9] + Reading frames:")
    for frame in gen_reading_frames(DNAStr):
        print("\t", frame)

    print("\n[10] + All proteins in 6 open reading frames:")
    for prot in all_proteins_from_orfs(DNAStr, 0, 0, True):
        print("\t", prot)

if __name__ == "__main__":
    main()