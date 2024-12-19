import os

# Listar directorios que empiezan con 'Conf-' y filtrar por los que tengan número mayor a 4
dirs = [i for i in os.listdir() if i.startswith('Conf-')]
dirs = [i for i in dirs if int(i.split('-')[1]) > 4]

# Iterar sobre los directorios
for dir in dirs:
    # Obtener archivo .gro y eliminar el salto de línea
    file_gro = os.popen(f'ls {dir}/*.gro').read().strip()
    
    # Verificar si el archivo .gro existe antes de proceder
    if not os.path.exists(file_gro):
        print(f"Archivo .gro no encontrado en {dir}")
        continue
    
    # Verificar si el archivo .mdp existe
    mdp_file = f"{dir}/MD_U.mdp"
    if not os.path.exists(mdp_file):
        print(f"Archivo .mdp no encontrado en {mdp_file}")
        continue
    
    # Verificar si el archivo .top existe
    top_file = f"{dir}/topol.top"
    if not os.path.exists(top_file):
        print(f"Archivo .top no encontrado en {top_file}")
        continue
    
    # Ejecutar el comando gmx_mpi grompp con los archivos correctos
    os.system(f'gmx_mpi grompp -c {file_gro} -f {mdp_file} -p {top_file} -o {dir}/{dir}_PMF.tpr -r {file_gro}')

