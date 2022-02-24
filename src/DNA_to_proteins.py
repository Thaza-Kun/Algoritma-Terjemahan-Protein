from typing import List
from structures import *

def validate_sequence(dna_seq: str) -> str:
    """Semak sama ada input adalah jujukan DNA adalah sah

    Args:
        dna_seq (str): Jujukan input

    Raises:
        ValueError: Akan muncul kalau ada huruf yang tidak dikenali dalam jujukan

    Returns:
        str: Jujukan yang sah
    """
    check_seq = dna_seq.upper()
    for nuc in check_seq:
        if nuc not in DNA_Bases:
            raise ValueError(f"Jujukan tidak sah. `{nuc}` bukan ahli {DNA_Bases=}")
    return check_seq

def make_complement(dna_seq: str) -> str:
    """Menghasilkan pasangan jujukan yang melekat bersama dengan jujukan input

    Args:
        dna_seq (str): Jujukan input

    Returns:
        str: Jujukan pasangan
    """
    mapping = str.maketrans('ATGC', 'TACG')
    return dna_seq.translate(mapping)

def make_reverse_complement(dna_seq: str) -> str:
    """Menterbalikkan pasangan jujukan (sebab pasangannya dibaca berlawanan arah)

    Args:
        dna_seq (str): Jujukan input

    Returns:
        str: Jujukan pasangan terbalik
    """
    return make_complement(dna_seq)[::-1]

def make_translation(dna_seq: str, init_pos: int = 0) -> list:
    """Menterjemahkan jujukan input menjadi senarai asid amino

    Args:
        dna_seq (str): Jujukan input
        init_pos (int, optional): Kedudukan bes pertama yang dibaca. Nilai asalnya 0.

    Returns:
        list: Senarai asid amino
    """
    acids = list()
    # Lelaran akan melangkau setiap tiga huruf sebab itulah panjang satu kodon
    for pos in range(init_pos, len(dna_seq) - 2, 3):
        codon = dna_seq[pos:pos + 3]
        acids.append(DNA_Codons_Mapping[codon])
    return acids

def make_reading_frames(dna_seq: str) -> List[List[str]]:
    """Menghasilkan kesemua enam (6) rangka bacaan daripada jujukan bes DNA yang diberikan.

    Args:
        dna_seq (str): Jujukan input

    Returns:
        List[List[str]]: Senarai rangka bacaan yang terdiri daripada senarai asid amino
    """
    frames = []
    frames.append(make_translation(dna_seq, 0))
    frames.append(make_translation(dna_seq, 1))
    frames.append(make_translation(dna_seq, 2))
    frames.append(make_translation(make_reverse_complement(dna_seq), 0))
    frames.append(make_translation(make_reverse_complement(dna_seq), 1))
    frames.append(make_translation(make_reverse_complement(dna_seq), 2))
    return frames

def make_proteins_from_one_reading_frame(amino_acid_seq: list) -> list:
    """Menghasilkan calon protein dengan mengenal pasti kodon pemula dan kodon berhenti

    Args:
        amino_acid_seq (list): Senarai asid amino

    Returns:
        list: Senarai calon protein
    """
    current_protein = []
    proteins = []
    for amino_acid in amino_acid_seq:
        if amino_acid == "_":
            # Hentikan rantain
            if current_protein:
                # Dapatkan rantaian jika wujud
                for p in current_protein:
                    proteins.append(p)
                current_protein = []
        else:
            if amino_acid == "M":
                # Mulakan rantaian
                current_protein.append("") 
            for i in range(len(current_protein)):
                # Operasi ini tidak akan buat apa-apa kalau tiada `""` yang dihasilkan oleh syarat sebelum ini
                current_protein[i] += amino_acid
    return proteins

def make_proteins_from_all_reading_frames(dna_seq: str, ordered=False) -> List[str]:
    """Menghasilkan calon protein daripada semua rangka bacaan yang dihasilkan dari sejujuk DNA

    Args:
        dna_seq (str): Jujukan input
        ordered (bool, optional): Sama ada mahu susun atau tak. Nilai asal False.

    Returns:
        List[str]: Senarai calon protein
    """
    frames = make_reading_frames(dna_seq)

    potential_proteins = []
    for frame in frames:
        proteins = make_proteins_from_one_reading_frame(frame)
        for p in proteins:
            potential_proteins.append(p)

    if ordered:
        potential_proteins = sorted(potential_proteins, key=len, reverse=True)
        
    return potential_proteins