import os
import subprocess

def convert_tga_to_png(directory="."):
    """
    Converts all .tga files in the specified directory to .png format.
    It checks if a .png file with the same name already exists before converting.
    """
    if not os.path.isdir(directory):
        print(f"Error: Directory '{directory}' not found.")
        return

    print(f"Searching for .tga files in '{os.path.abspath(directory)}'...")
    converted_count = 0

    for filename in os.listdir(directory):
        if filename.lower().endswith(".tga"):
            base_name = os.path.splitext(filename)[0]
            tga_path = os.path.join(directory, filename)
            png_filename = base_name + ".png"
            png_path = os.path.join(directory, png_filename)

            if os.path.exists(png_path):
                print(f"Skipping '{filename}': '{png_filename}' already exists.")
            else:
                print(f"Converting '{filename}' to '{png_filename}'...")
                try:
                    # Use 'convert' command line tool
                    subprocess.run(["convert", tga_path, png_path], check=True)
                    print(f"Successfully converted '{filename}'.")
                    converted_count += 1
                except subprocess.CalledProcessError as e:
                    print(f"Error converting '{filename}': {e}")
                except FileNotFoundError:
                    print("Error: 'convert' command not found. Make sure ImageMagick is installed and in your system's PATH.")
                    print("You can usually install it via your package manager (e.g., 'sudo apt-get install imagemagick' on Ubuntu, 'brew install imagemagick' on macOS).")
                    return

    if converted_count == 0:
        print("No new .tga files were found or converted.")
    else:
        print(f"\nConversion complete! Converted {converted_count} .tga files to .png.")

# --- How to use the script ---
if __name__ == "__main__":
    # You can specify the directory where your .tga files are.
    # If you leave it blank or use '.', it will use the current directory where the script is run.
    target_directory = "." # Change this to your directory if needed, e.g., "C:/Users/YourName/MyRenders"
    convert_tga_to_png(target_directory)
