import os

# Definimos el tamaño límite en bytes (1 GB)
SIZE_LIMIT = 1 * 1024 * 1024 * 1024  # 1 GB

# Directorio base donde están las carpetas Conf-*
base_directory = os.getcwd()

# Recorremos todas las carpetas en el directorio base
for folder in os.listdir(base_directory):
    if folder.startswith("Conf-") and os.path.isdir(folder):
        folder_path = os.path.join(base_directory, folder)
        
        # Buscamos un archivo con extensión .xtc
        xtc_files = [f for f in os.listdir(folder_path) if f.endswith(".xtc")]

        if xtc_files:
            # Si hay más de un archivo .xtc, tomamos el primero (suponemos uno por carpeta)
            xtc_file = xtc_files[0]
            xtc_path = os.path.join(folder_path, xtc_file)

            # Obtenemos el tamaño del archivo .xtc
            xtc_size = os.path.getsize(xtc_path)

            if xtc_size > SIZE_LIMIT:
                print(f"El archivo {xtc_file} en {folder} supera 1GB ({xtc_size / (1024 * 1024):.2f} MB).")

                # Buscamos un archivo con extensión .trr
                trr_files = [f for f in os.listdir(folder_path) if f.endswith(".trr")]

                if trr_files:
                    # Si hay más de un archivo .trr, tomamos el primero (suponemos uno por carpeta)
                    trr_file = trr_files[0]
                    trr_path = os.path.join(folder_path, trr_file)

                    # Eliminamos el archivo .trr
                    try:
                        os.remove(trr_path)
                        print(f"Archivo {trr_file} eliminado en {folder}.")
                    except Exception as e:
                        print(f"Error al eliminar {trr_file} en {folder}: {e}")
                else:
                    print(f"No se encontró ningún archivo .trr en {folder}.")
            else:
                print(f"El archivo {xtc_file} en {folder} no supera 1GB ({xtc_size / (1024 * 1024):.2f} MB). No se eliminará nada.")
        else:
            print(f"No se encontró ningún archivo .xtc en {folder}.")

