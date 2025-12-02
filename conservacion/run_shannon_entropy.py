import numpy as np
import argparse
import os


iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"] 

# dictionary to map from amino acid to its row/column in a similarity matrix
aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']

# dictionary to map from amino acid to its row/column in a similarity matrix
aa_to_index = {aa: i for i, aa in enumerate(aas)}

def calculate_sequence_weights(msa):
    """ 
    Devuelve una lista de pesos por cada secuencia en un MSA 
    según el método de Henikoff '94, donde las secuencias 
    más diferentes tienen un mayor peso. 
    """

    seq_weights = [0. for _ in range(len(msa))]
    
    for i in range(len(msa[0])):
        freq_counts = [0 for _ in range(len(aas))]
        
        for j in range(len(msa)):
            if msa[j][i] in aa_to_index:  # Verifica si el aminoácido está en el diccionario
                freq_counts[aa_to_index[msa[j][i]]] += 1

        num_observed_types = sum(1 for count in freq_counts if count > 0)

        for j in range(len(msa)):
            if msa[j][i] in aa_to_index:  # Verifica si el aminoácido está en el diccionario
                d = freq_counts[aa_to_index[msa[j][i]]] * num_observed_types
                if d > 0:
                    seq_weights[j] += 1. / d

    for w in range(len(seq_weights)):
        seq_weights[w] /= len(msa[0])

    return seq_weights

def read_fasta_alignment(filename):
    """
    Lee el MSA y devuelve dos listas;
    identifiers y secuencia.
    """

    iupac_alphabet = "ACDEFGHIKLMNPQRSTVWY"
    names = []
    alignment = []

    with open(filename) as f:
        cur_seq = ''
        for line in f:
            line = line.strip()
            if line.startswith(';') or len(line) == 0:
                continue
            if line.startswith('>'):
                if cur_seq: #me fijo si tengo algo escrito en mi secuencia para procesarlo
                    cleaned_seq = ''.join([aa if aa in iupac_alphabet else '-' for aa in cur_seq])
                    alignment.append(cleaned_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
                names.append(line[1:]) #apendeo el nombre de la linea con '>' 
                cur_seq = ''
            else:
                cur_seq += line

        if cur_seq: #sirve para la ultima secuencia
            cleaned_seq = ''.join([aa if aa in iupac_alphabet else '-' for aa in cur_seq])
            alignment.append(cleaned_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))

    return names, alignment
import math

def weighted_gap_penalty(col, seq_weights):
    """ Calculate the simple gap penalty multiplier for the column. If the 
    sequences are weighted, the gaps, when penalized, are weighted 
    accordingly. """

    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1.] * len(col)
    
    gap_sum = 0.
    for i in range(len(col)):
        if col[i] == '-':
            gap_sum += seq_weights[i]

    return 1 - (gap_sum / sum(seq_weights))

def gap_percentage(col):
    """Return the percentage of gaps in col."""
    num_gaps = 0.

    for aa in col:
        if aa == '-': num_gaps += 1

    return num_gaps / len(col)


def frecuencia_pesada(col, seq_weights, pc_amount):
    """ me devuelve la frecuencia de los aa multiplicadas por el peso
    segun la secuencia que correspondan usando pseudocount para los AA
    que no tengan frecuencia."""
    # si no hay matcheo uso el mismo peso para todas las secuencias
    if len(seq_weights) != len(col):
        seq_weights = [1. for i in range(len(col))]

    aa_num = 0
    freq_counts = [pc_amount for i in range(len(amino_acids))] # respeta el orden de lso aminoacidos

    for aa in amino_acids:
        for j in range(len(col)):
            if col[j] == aa:
                freq_counts[aa_num] += 1 * seq_weights[j]

        aa_num += 1

    for j in range(len(freq_counts)):
        freq_counts[j] = freq_counts[j] / (sum(seq_weights) + len(amino_acids) * pc_amount)

    return freq_counts

def get_column(col_num, alignment):
    """
    Devuelve una columna del MSA como lista
    """
    col = []
    for seq in alignment:
        if col_num < len(seq):
            col.append(seq[col_num])
    return col
def window_score(scores, window_len, lam=.5):
    """
    This function takes a list of scores and a length and transforms them 
    so that each position is a weighted average of the surrounding positions. 
    Positions with scores less than zero are not changed and are ignored in the 
    calculation. Here window_len is interpreted to mean window_len residues on 
    either side of the current residue
    """

    w_scores = scores[:]

    for i in range(window_len, len(scores) - window_len):
        if scores[i] < 0: 
            continue

        sum = 0.
        num_terms = 0.
        for j in range(i - window_len, i + window_len + 1):
            if i != j and scores[j] >= 0:
                num_terms += 1
                sum += scores[j]

        if num_terms > 0:
            w_scores[i] = (1 - lam) * (sum / num_terms) + lam * scores[i]

    return w_scores
amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"] 
blosum_background_distr = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]
bg_distribution = blosum_background_distr[:]


def conservation(msa, seq_pesos, funcion_conservacion):
    scores = []
    for i in range(len(msa[0])):
        columna = get_column(col_num=i, alignment=msa)
        a = funcion_conservacion(columna,bg_distribution, seq_pesos) #js_divergence(columna, bg_distribution, pesos)
        scores.append(a)
    return(scores)


def obtener_indices_no_guion(cadena):
    indices = [i for i, char in enumerate(cadena) if char != '-']
    return indices


### funciones de conservacion
def js_divergence(col, bg_distr, seq_weights, gap_penalty=1):
    """
    Me de vuelde el valore de Divergencia de Jensen-Shannon Divergence para la columna usando 
    los alineamientos de BLOSUM62 como frecuencia de fondos.
    """
    distr = bg_distr[:]
    PSEUDOCOUNT= 1e-6


    fc = frecuencia_pesada(col, seq_weights, PSEUDOCOUNT)

    if len(fc) != len(distr): return -1

    #  hago la distribucion r
    r = []
    for i in range(len(fc)):
        r.append(.5 * fc[i] + .5 * distr[i])

    d = 0.
    for i in range(len(fc)):
        if r[i] != 0.0:
            if fc[i] == 0.0:
                d += distr[i] * math.log(distr[i]/r[i], 2)
            elif distr[i] == 0.0:
                d += fc[i] * math.log(fc[i]/r[i], 2) 
            else:
                d += fc[i] * math.log(fc[i]/r[i], 2) + distr[i] * math.log(distr[i]/r[i], 2)

    # d /= 2 * math.log(len(fc))
    d /= 2

    if gap_penalty == 1: 
        return d * weighted_gap_penalty(col, seq_weights)
    else: 
        return d

def property_relative_entropy(col, bg_distr, seq_weights, gap_penalty=1):
    """Calculate the relative entropy of a column col relative to a
    partition of the amino acids. Similar to Williamson '95. sim_matrix is
    ignored, but could be used to define the sets. See shannon_entropy()
    for more general info. """
    PSEUDOCOUNT= 1e-6
    # Mirny and Shakn. '99
    #property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]

    # Williamson '95
    property_partition = [['V','L', 'I','M'], ['F','W','Y'], ['S','T'], ['N','Q'], ['H','K','R'], ['D','E'], ['A','G'], ['P'], ['C']] #LE SAQUE LOS GAPS MIRAR QUE ONDA SI FRECUENCIA PESADA LO DEVUELVE COMO FREC LA CANTIDAD DE GAPS

    prop_bg_freq = []
    if len(bg_distr) == len(property_partition):
        prop_bg_freq = bg_distr
    else: prop_bg_freq = [0.248, 0.092, 0.114, 0.075, 0.132, 0.111, 0.161, 0.043, 0.024, 0.000]
    fc = frecuencia_pesada(col, seq_weights, PSEUDOCOUNT) # sum the aa frequencies to get the property frequencies
    prop_fc = [0.] * len(property_partition)
    for p in range(len(property_partition)):
        for aa in property_partition[p]:
            prop_fc[p] += fc[aa_to_index[aa]]

    d = 0. 
    for i in range(len(prop_fc)):
        if prop_fc[i] != 0 and prop_bg_freq[i] != 0:d += prop_fc[i] * math.log(prop_fc[i] / prop_bg_freq[i], 2)
    if gap_penalty == 1: return d * weighted_gap_penalty(col, seq_weights)
    else:
        return(d)

def shannon_entropy(col, seq_weights, gap_penalty=1):
    """Calculates the Shannon entropy of the column col. sim_matrix  and 
    bg_distr are ignored. If gap_penalty == 1, then gaps are penalized. The 
    entropy will be between zero and one because of its base. See p.13 of 
    Valdar 02 for details. The information score 1 - h is returned for the sake 
    of consistency with other scores."""
    PSEUDOCOUNT= 1e-6
    fc = frecuencia_pesada(col, seq_weights, PSEUDOCOUNT)

    h = 0. 
    for i in range(len(fc)):
        if fc[i] != 0:
            h += fc[i] * math.log(fc[i])

#    h /= math.log(len(fc))
    h /= math.log(min(len(fc), len(col)))

    inf_score = 1 - (-1 * h)
    if gap_penalty == 1: return inf_score * weighted_gap_penalty(col, seq_weights)
    else: 
	    return (inf_score)
def property_entropy(col, seq_weights, gap_penalty=1):
    """Calculate the entropy of a column col relative to a partition of the 
    amino acids. Similar to Mirny '99. sim_matrix and bg_distr are ignored, but 
    could be used to define the sets. """

    # Mirny and Shakn. '99
    #property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]

    # Williamson '95
    property_partition = [['V','L', 'I','M'],['F','W','Y'], ['S','T'], ['N','Q'], ['H','K','R'], ['D','E'], ['A','G'], ['P'], ['C']]
    PSEUDOCOUNT= 1e-6
    fc = frecuencia_pesada(col, seq_weights, PSEUDOCOUNT)

    # sum the aa frequencies to get the property frequencies
    prop_fc = [0.] * len(property_partition)
    for p in range(len(property_partition)):
        for aa in property_partition[p]: prop_fc[p] += fc[aa_to_index[aa]]
    h = 0.0 
    for i in range(len(prop_fc)):
        if prop_fc[i] != 0: h += prop_fc[i] * math.log(prop_fc[i])
    h /= math.log(min(len(property_partition), len(col)))

    inf_score = 1 - (-1 * h)

    if gap_penalty == 1: return inf_score * weighted_gap_penalty(col, seq_weights)
    else: 
	    return inf_score




def run_conservation_shannon_entropy(path_msa):
    print('La Primer Seq del MSA tiene que ser nuestra referencia')
    names_msa, msa = read_fasta_alignment(path_msa)
    pesos = calculate_sequence_weights(msa)
    scores = conservation(msa, pesos, shannon_entropy)
    indices_columnas = obtener_indices_no_guion(msa[0])
    wc_df = []
    for i in indices_columnas:
        #wc_df.append(wc_conservation[i])
        numero_formateado = float("{:.2f}".format(scores[i]))
        wc_df.append(numero_formateado)
    return wc_df



if __name__ == '__main__':
    # Configuración de los argumentos de línea de comandos
    parser = argparse.ArgumentParser(
        description="Calcula los scores de conservación (Entropía de Shannon) para un alineamiento FASTA."
    )
    
    # Argumento para el archivo de entrada
    parser.add_argument(
        'input_file', 
        type=str, 
        help='Ruta al archivo de alineamiento FASTA (.in)'
    )
    
    # Argumento para el archivo de salida
    parser.add_argument(
        'output_file', 
        type=str, 
        help='Ruta al archivo donde se guardarán los scores (.txt)'
    )
    
    args = parser.parse_args()
    
    try:
        # 1. Ejecutar la lógica de conservación
        scores_result = run_conservation_shannon_entropy(args.input_file)
        
        # 2. Guardar el resultado en el archivo de salida
        np.savetxt(args.output_file, scores_result, fmt='%.2f')
        
        print('--------------------------------------------------')
        print(f"Éxito: Scores de conservación guardados en: {args.output_file}")
        
    except FileNotFoundError:
        print(f"ERROR: No se encontró el archivo de entrada: {args.input_file}")
    except Exception as e:
        print(f"Ocurrió un error inesperado durante la ejecución: {e}")
