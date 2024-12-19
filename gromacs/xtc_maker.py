import os
import subprocess

def run_gmx_command():
    # Obtener el directorio actual
    current_dir = os.getcwd()

    # Iterar por las carpetas que coincidan con el formato Conf-{número}
    for folder in os.listdir(current_dir):
        if os.path.isdir(folder) and folder.startswith("Conf-"):
            folder_path = os.path.join(current_dir, folder)

            # Construir los nombres de los archivos basados en la carpeta actual
            trr_file = os.path.join(folder_path, f"{folder}_PMF.trr")
            tpr_file = os.path.join(folder_path, f"{folder}_PMF.tpr")
            output_file = os.path.join(folder_path, f"{folder}_PMF_noPBC.xtc")

            # Construir el comando
            command = f"echo 2 0 | gmx_mpi trjconv -f {trr_file} -s {tpr_file} -o {output_file} -pbc mol -center"

            print(f"Processing: {folder}")

            # Ejecutar el comando
            try:
                subprocess.run(command, shell=True, check=True)
                print(f"Successfully processed {folder}")
            except subprocess.CalledProcessError as e:
                print(f"Error processing {folder}: {e}")

if __name__ == "__main__":
    run_gmx_command()

