import json
import os 
import sys

comandos_ayuda = ['?','-h','-help','-H','h','help', '']
# Abre el archivo JSON en modo lectura


os.system('')
def fasta_maker_from_hits_FSeek(archivo_json):

    with open(archivo_json, 'r') as DATA:
        datos = json.load(DATA)
    #QUERY
    header_reference = 'reference'
    seq_reference = datos[0]['query']['sequence']
    with open(f'homologues_foldseek.fasta', 'w') as w:
        w.write(f'>{header_reference}\n{seq_reference}\n')
        for linea in (datos[0]['results'][0]['alignments']):
            inicio = linea['dbStartPos']; fin = linea['dbEndPos']
            seq = linea['tSeq'] #guardo la secuencia del hit
            seq_hit = seq[inicio-1:fin] #guardo solo la parte que alinea contra el query
            header = linea['target']

            w.write(f'>{header}\n{seq_hit}\n')
    w.close()

        
if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] in comandos_ayuda:
        print('Uso del programa: \npython3 programa.py archivo.json')
        sys.exit(1)
    archivo_json = sys.argv[1]
    fasta_maker_from_hits_FSeek(archivo_json)  
