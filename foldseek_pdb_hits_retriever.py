#Este archivo toma el .json de salida de foldseek search y retribuye los hits, considerando
#solo las cadenas que alinearon y tambien trae ligandos asociados
import json
import os 
import sys
import subprocess
import urllib
from Bio.PDB import PDBParser, PDBIO, Select
import tqdm
comandos_ayuda = ['?','-h','-help','-H','h','help', '']
# Abre el archivo JSON en modo lectura

comandos_ayuda = ['?','-h','-help','-H','h','help', '']
# Abre el archivo JSON en modo lectura


def load_json(archivo_json):

    with open(archivo_json, 'r') as DATA:
        datos = json.load(DATA)
    return datos

def GetPDB(pdbid):

    downloadurl="https://files.rcsb.org/download/"
    file = pdbid+'.pdb'
    url = downloadurl + file

    command = ["wget", url, "-O", file]
    subprocess.run(command)
  
class Selector(Select):

    def __init__(self, cadena):
        self.chain = cadena[0]

    def accept_residue(self, res):
        statement = (res.get_full_id()[3][0].startswith('H')) or (res.get_full_id()[2] == self.chain)
        return statement

def pdb_hit_fseek_retriever(datos):
    datos = load_json(datos)

    for linea in (datos[0]['results'][0]['alignments']):
        pdbid, chain = (linea['target'].split('_'))[0], (linea['target'].split('_'))[1]
        try:
            GetPDB(pdbid)
            parser = PDBParser(QUIET=True)
            directorio = subprocess.getoutput('pwd')
            file_pdb = directorio+'/'+pdbid+'.pdb'
            structure = parser.get_structure('estructura', file_pdb)
            selector = Selector('A')
            io = PDBIO()
            io.set_structure(structure)
            io.save(pdbid+'_curated.pdb', selector)
        except: print('no encontre a ', pdbid)



if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] in comandos_ayuda:
        print('Uso del programa: \npython3 programa.py archivo.json')
        sys.exit(1)
    archivo_json = sys.argv[1]
    pdb_hit_fseek_retriever(archivo_json)