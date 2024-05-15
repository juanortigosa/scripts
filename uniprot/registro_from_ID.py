from Bio import SeqIO
from Bio import ExPASy
def uniprot_handler(uniprotID): 
  '''
    con el UniProtID elegimos la secuencia ----> Ej. UGGT 1 humana = Q9NYU2, myoglobin human = P02144 
    '''
  with ExPASy.get_sprot_raw(uniprotID) as handle:
    record = SeqIO.read(handle, 'swiss')
    return record
#  name = record.description.split(';')[0].split('=')[1]
#  seq_query = record.seq
#  name_query = record.description.split(';')[0].split('=')[1]
#  #creo un archivo fasta para hacer un BLAST contra swissprot/uniprotKB
#  with open(f'{archivo_salida_fasta}', 'w') as writable:
#    writable.write('>'+name_query+'\n'+str(seq_query))
#    writable.close()
#  return(name_query, seq_query)