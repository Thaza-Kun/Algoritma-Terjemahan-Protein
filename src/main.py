from os import path, getcwd

from DNA_to_proteins import *
from utilities import color_text_on_print, read_FASTA

DATA_PATH = 'data'

def main(file: str):
    """Baca data dari fail kemudian transkripsikan semua kemungkinan protein

    Args:
        file (str): fail yang mengandungi data
    """
    get_extension = file.split('.')[-1]
    if get_extension == 'fasta':
        file = path.join(getcwd(), DATA_PATH, file)
        FASTAData = read_FASTA(file)
        DNAList = []
        for i in list(FASTAData.values()):
            DNAList.append(validate_sequence(i))
        DNAStr = DNAList[0]
    else:
        file = path.join(getcwd(), DATA_PATH, file)
        with open(file, 'r') as f:
            DNAStr = f.read()

    print(f"\nJujukan input ({len(DNAStr)} bes):\n{color_text_on_print(DNAStr)}\n")

    print("Semua senarai asid amino:")
    for num, frame in enumerate(make_reading_frames(DNAStr)):
        print(f"{num + 1}.", ''.join(frame))

    print("\nSemua protein yang mungkin dihasilkan:")
    for prot in make_proteins_from_all_reading_frames(DNAStr, 0, 0, True):
        print("-", prot)

if __name__ == "__main__":
    main('homo-sapiens-ins-var2.fasta')