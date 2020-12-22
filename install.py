import subprocess
import argparse
import os.path
import platform

def extract_all_zips_in_folder(folder):
    print("extracting files in " + folder)
    if not folder.endswith("/"):
        folder += "/"

    # extract all zip files
    with os.scandir(folder) as files:
        escaped_folder = folder.replace('"', '\\"')
        for entry in files:
            if entry.name.endswith(".zip"):
                print("Extracting " + entry.name)
                subprocess.call('unzip -n -d "' + escaped_folder + '" "' + escaped_folder + entry.name.replace('"', '\\"') + '"', shell=True)


print("\n**************************************************")
print("                       COM2POSE                     ")
print("                  Installation Script               ")
print("**************************************************")

subprocess.check_call(['Rscript', 'ext/TEPIC/TEPIC/Code/installRpackages.R'], shell=False)
