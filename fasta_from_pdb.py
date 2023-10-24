from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

codigo_aminoacidos = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V"
}

def pdb_to_fasta(pdb_file, fasta_file):
    # Leer el archivo PDB y extraer la secuencia de aminoácidos
    parser = PDBParser(QUIET=True)
    estructura = parser.get_structure("estructura", pdb_file)
    model = estructura[0]
    sequence = ""

    for chain in model:
    	for residue in chain:
    		sequence += codigo_aminoacidos[residue.get_resname()]
  
    # Crear un registro de secuencia FASTA
    record = SeqRecord(Seq(sequence), id=f"{fasta_file}")

    # Escribir el registro de secuencia FASTA en un archivo
    with open(fasta_file, 'w') as output_file:
        SeqIO.write(record, output_file, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python programa.py archivo.pdb salida.fasta")
        sys.exit(1)
    pdb_file = sys.argv[1]
    fasta_file = sys.argv[2]
    pdb_to_fasta(pdb_file, fasta_file)

