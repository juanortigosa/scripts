import argparse
import os

# Escala de Hidrofobicidad de Kyte & Doolittle (1982)
# Valores de -4.5 (más hidrofílico) a +4.5 (más hidrofóbico)
# Fuente: https://web.expasy.org/protscale/
HYDROPHOBICITY_SCALE = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLU': -3.5, 'GLN': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}

def assign_hydrophobicity_to_bfactors(input_pdb_path: str, output_pdb_path: str):
    """
    Lee un archivo PDB, reemplaza la columna de B-factor (temperatura factor)
    con el valor de hidrofobicidad de Kyte & Doolittle del residuo
    correspondiente, y guarda el resultado en un nuevo archivo PDB.
    
    Args:
        input_pdb_path (str): Ruta al archivo PDB de entrada.
        output_pdb_path (str): Ruta al archivo PDB de salida.
    """
    
    if not os.path.exists(input_pdb_path):
        raise FileNotFoundError(f"Error: No se encontró el archivo de entrada en la ruta: {input_pdb_path}")

    # Diccionario para almacenar el valor de hidrofobicidad por ID de Residuo
    # Llave: (residuo_nombre, numero_residuo, cadena) -> Valor: hidrofobicidad
    residue_hydrophobicity = {}
    
    print(f"DEBUG: Procesando PDB de entrada: {os.path.basename(input_pdb_path)}")
    
    # 1. Primera pasada: Leer y asignar valores de hidrofobicidad a cada residuo único
    with open(input_pdb_path, 'r') as infile:
        for line in infile:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Campos relevantes en el formato PDB (basado en índice de columna):
                # Residue Name: columnas 17-20 (ej. 'ALA', 'GLY')
                # Chain ID: columnas 21 (ej. 'A', 'B')
                # Residue Number: columnas 22-26 (ej. ' 10 ', ' 123A')
                
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_num = line[22:26].strip()
                
                # Crear una clave única para cada residuo
                res_key = (res_name, res_num, chain_id)
                
                if res_key not in residue_hydrophobicity:
                    # Buscar la hidrofobicidad. Si no está, usar 0.0 y notificar.
                    try:
                        hydro_value = HYDROPHOBICITY_SCALE[res_name]
                        residue_hydrophobicity[res_key] = hydro_value
                    except KeyError:
                        # Esto maneja residuos no estándar (ej. modificados, HETATM si no es agua)
                        hydro_value = 0.0 # Valor neutro por defecto
                        residue_hydrophobicity[res_key] = hydro_value
                        # Se podría agregar una advertencia, pero por ahora se omite para limpieza
    
    print(f"DEBUG: Se asignaron valores de hidrofobicidad a {len(residue_hydrophobicity)} residuos únicos.")

    # 2. Segunda pasada: Escribir el nuevo archivo PDB con los B-factors actualizados
    with open(input_pdb_path, 'r') as infile, open(output_pdb_path, 'w') as outfile:
        
        # Mapear el rango de hidrofobicidad al rango deseado (0.00 a 1.00)
        all_values = [v for v in residue_hydrophobicity.values() if v != 0.0]
        if all_values:
            min_hydro = min(all_values)
            max_hydro = max(all_values)
            hydro_range = max_hydro - min_hydro
            
            def map_to_bfactor_range(value):
                if value == 0.0 or hydro_range == 0: 
                    return 0.00
                # Normaliza el valor entre 0 y 1.
                normalized = (value - min_hydro) / hydro_range
                return normalized * 1.00 # Multiplicar por 1.00 para mantener el rango 0-1
        else:
             def map_to_bfactor_range(value): return 0.00 # Fallback si solo hay un tipo de residuo
        
        # Reescritura línea por línea
        for line in infile:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_num = line[22:26].strip()
                res_key = (res_name, res_num, chain_id)
                
                # Obtener la hidrofobicidad asignada (valor original)
                hydro_value_orig = residue_hydrophobicity.get(res_key, 0.0)
                
                # Mapear el valor para que esté en el rango 0.00 - 1.00
                bfactor_value = map_to_bfactor_range(hydro_value_orig)
                
                # Formatear el B-factor (6 columnas, precisión de 2 decimales, F6.2).
                # Columna 61-66 (Python index 60-66).
                # Usamos rjust(6) para garantizar que el valor se alinee a la derecha.
                # Nota: Si bfactor_value es 1.00, la cadena es "  1.00". Si es 0.10, es "  0.10".
                bfactor_str = f"{bfactor_value:6.2f}"
                
                # Construir la nueva línea
                # line[:60] -> Todo hasta la columna 60 (incluye Occupancy)
                # bfactor_str -> Las 6 columnas del nuevo B-factor (61-66)
                # line[66:] -> El resto de la línea (Element, Charge, etc.)
                new_line = line[:60] + bfactor_str + line[66:]
                
                outfile.write(new_line)
            else:
                # Escribir todas las otras líneas (HEADER, REMARK, TER, etc.) sin modificar
                outfile.write(line)

    print(f"Éxito: El archivo modificado ha sido guardado en: {output_pdb_path}")
    print("--------------------------------------------------")
    print("Nota: El B-factor ahora está escalado al rango [0.00, 1.00] para su visualización.")


if __name__ == '__main__':
    # Configuración de los argumentos de línea de comandos
    parser = argparse.ArgumentParser(
        description="Reemplaza los B-factors de un archivo PDB por los valores de hidrofobicidad de Kyte & Doolittle.",
        epilog="Ejemplo de uso: python pdb_hydrophobicity.py 1abc.pdb 1abc_hydro.pdb"
    )
    
    # Argumento para el archivo de entrada
    parser.add_argument(
        'input_pdb', 
        type=str, 
        help='Ruta al archivo PDB de entrada (ej. input.pdb).'
    )
    
    # Argumento para el archivo de salida
    parser.add_argument(
        'output_pdb', 
        type=str, 
        help='Ruta al archivo PDB de salida donde se guardará el resultado (ej. output.pdb).'
    )
    
    args = parser.parse_args()
    
    try:
        assign_hydrophobicity_to_bfactors(args.input_pdb, args.output_pdb)
        
    except FileNotFoundError as e:
        print(e)
    except Exception as e:
        print(f"Ocurrió un error inesperado: {e}")
