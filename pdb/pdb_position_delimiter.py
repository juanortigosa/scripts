# Delimitador de posiciones
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
    import sys

    if len(sys.argv) != 5:
        if sys.argv[1] in ["-h", "--help", "-h", "--help"]:
            print("Uso: python programa.py archivo.pdb inicio fin archivo_salida",
            "Este programa recorta un pdb segun posiciones")
            sys.exit(0)
        else:
            print("Uso incorrecto. Para ayuda, ejecute con -h, --help, -h, o --help.")
            sys.exit(1)

    archivo_entrada = sys.argv[1]
    residuo_menor = int(sys.argv[2])
    residuo_mayor = int(sys.argv[3])
    archivo_salida = sys.argv[4]

    seleccionar_y_guardar_rango_pdb(archivo_entrada, residuo_menor, residuo_mayor, archivo_salida)

if __name__ == "__main__":
    main()

