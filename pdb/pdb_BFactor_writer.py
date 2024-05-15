import argparse
import numpy as np
import sys
from Bio.PDB import PDBIO, PDBParser

import argparse
import numpy as np
import sys
from Bio.PDB import PDBIO, PDBParser

def make_file(pdb_file, b_factors_file, output_file):
    b_factors_file = np.loadtxt(b_factors_file)
    file = ''
    # Lee la primera línea para obtener el primer valor de B-score
    pdb_file = open(pdb_file).readlines()
    numero_antiguo = None
    indice_bscore = 0
    for line in pdb_file:
        if line.startswith('ATOM'):
            numero = int(line[22:26])
            atomo = line.split()[2][0]
            if numero == numero_antiguo or numero_antiguo == None:
                # Si es el mismo residuo, asigna el mismo B-score
                line = line[:60] + str("{:6.2f}".format(100 * b_factors_file[indice_bscore])) + line[66:77] + atomo + '   \n'
            else:
                indice_bscore += 1
                line = line[:60] + str("{:6.2f}".format(100 * b_factors_file[indice_bscore])) + line[66:77] + atomo + '   \n'
            numero_antiguo = numero
            file += line
    # Escribir el registro de secuencia FASTA en un archivo
    with open(output_file, 'w') as output:
        output.write(file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python programa.py archivo.pdb BFactors(por-residuo).txt output.pdb")
        sys.exit(1)
    pdb_file = sys.argv[1]
    b_factors_file = sys.argv[2]
    output_file = sys.argv[3]  # Añade esta línea para definir output_file
    make_file(pdb_file, b_factors_file, output_file)
