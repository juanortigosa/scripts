import os
import subprocess

def run_exceptions(base_dir):
    # Lista todas las carpetas en el directorio base
    folders = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))]

    # Filtrar solo carpetas que cumplan con "Conf-17" o superiores
    for folder in folders:
        if folder.startswith("Conf-"):
            try:
                number = int(folder.split("-")[1])
                if number >= 17:
                    folder_path = os.path.join(base_dir, folder)

                    # Nombres específicos para archivos en Conf-{number}
                    trr_file = os.path.join(folder_path, f"{folder}_PMF.trr")
                    tpr_file = os.path.join(folder_path, f"{folder}_PMF.tpr")
                    output_file = os.path.join(folder_path, f"{folder}_PMF_noPBC.xtc")

                    # Comando
                    command = (
                        f"echo 2 0 | gmx_mpi trjconv -f {trr_file} -s {tpr_file} "
                        f"-o {output_file} -pbc mol -center"
                    )

                    print(f"Processing: {folder}")
                    subprocess.run(command, shell=True, check=True)
            except ValueError:
                print(f"Skipping invalid folder name: {folder}")

if __name__ == "__main__":
    base_directory = os.getcwd()  # Cambiar si es necesario
    run_exceptions(base_directory)

