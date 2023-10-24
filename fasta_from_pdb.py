from Bio import PDB
from Bio.SeqRecord import SeqRecord
import sys

def pdb_to_fasta(pdb_file, fasta_file):
    # Leer el archivo PDB y extraer la secuencia de aminoácidos
    sequence = ""
    with open(pdb_file, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") and line[17:20] == "CA ":
                residue = line[21:23]
                sequence += residue

    # Crear un registro de secuencia FASTA
    record = SeqRecord(sequence, id="PDB_Seq")

    # Escribir el registro de secuencia FASTA en un archivo
    with open(fasta_file, 'w') as output_file:
        SeqIO.write(record, output_file, "fasta")

if __name__ == "__main__":

    if sys.argv[1] == "-h" or sys.argv[1] == '-help' or sys.argv[1] == '-h' or sys.argv[1] == '-H':
        print('Mira capo se usa asi esto => python programa.py archivo.pdb salida.fasta')
        
        
    elif len(sys.argv) != 3:
        print("Uso: python programa.py archivo.pdb salida.fasta")
        sys.exit(1)

    else: 
        pdb_file = sys.argv[1]
        fasta_file = sys.argv[2]

        pdb_to_fasta(pdb_file, fasta_file)
