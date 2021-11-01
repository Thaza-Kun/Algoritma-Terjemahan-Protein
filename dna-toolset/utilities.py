"""
Author  : Murthadza bin Aznam
Date    : 2021-11-01

This file contains helping functions.
"""

def colored(dna_seq: str):
    bcolors = {
        "A": "\033[92m",
        "C": "\033[94m",
        "G": "\033[93m",
        "T": "\033[91m",
        "U": "\033[91m",
        "reset": "\033[0;0m",
    }

    tempStr = ""

    for nuc in dna_seq:
        if nuc in bcolors:
            tempStr += bcolors[nuc] + nuc
        else:
            tempStr += bcolors["reset"] + nuc
    return tempStr + bcolors["reset"]

def readFASTA(filePath) -> dict:
    with open(filePath, 'r') as f:
        FASTAFile = [l.strip() for l in f.readlines()]
    FASTADict = {}
    FASTALabel = ""
    for line in FASTAFile:
        if ">" in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line
    return FASTADict

def test():
    print(readFASTA("./dataset-rosalind/test-fasta.txt"))

if __name__=="__main__":
    test()