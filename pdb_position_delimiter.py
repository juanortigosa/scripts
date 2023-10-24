#Delimitador de posiciones
import argparse
from Bio.PDB import PDBParser, PDBIO, Select

class RangoSelector(Select):
    def __init__(self, inicio, fin):
        self.inicio = inicio
        self.fin = fin

    def accept_residue(self, residue):
        numero_residuo = residue.id[1]
        return self.inicio <= numero_residuo <= self.fin

def seleccionar_y_guardar_rango_pdb(input_pdb, inicio, fin, output_pdb):
    parser = PDBParser(QUIET=True)
    estructura = parser.get_structure("estructura", input_pdb)

    # Crea un objeto Selector para el rango deseado
    selector = RangoSelector(inicio, fin)

    # Guarda la estructura seleccionada en un nuevo archivo PDB
    io = PDBIO()
    io.set_structure(estructura)
    io.save(output_pdb, selector)

def main():
    parser = argparse.ArgumentParser(description="Selecciona y guarda un rango de posiciones de un archivo PDB.")
    parser.add_argument("archivo_entrada", help="Archivo PDB de entrada")
    parser.add_argument("residuo_menor", type=int, help="Número del residuo menor en el rango")
    parser.add_argument("residuo_mayor", type=int, help="Número del residuo mayor en el rango")
    parser.add_argument("archivo_salida", help="Nombre del archivo de salida PDB")

    args = parser.parse_args()
    archivo_entrada = args.archivo_entrada
    residuo_menor = args.residuo_menor
    residuo_mayor = args.residuo_mayor
    archivo_salida = args.archivo_salida

    seleccionar_y_guardar_rango_pdb(archivo_entrada, residuo_menor, residuo_mayor, archivo_salida)

if __name__ == "__main__":
    main()
