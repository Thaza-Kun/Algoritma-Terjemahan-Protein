def color_text_on_print(dna_seq: str):
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

def read_FASTA(filePath: str) -> dict:
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