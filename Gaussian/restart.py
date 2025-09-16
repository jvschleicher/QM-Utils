import glob
import os
import shutil
from pathlib import Path

path = os.getcwd()

log_files = []

version = input("Which version is the original that needs to be checked? e.g. leave blank for first iteration, _v2, _v3, ... \n ")

version_new = input("Which version needs to be created? e.g. _v2, _v3, ... \n")

end = version + ".log"

for file in os.listdir(path):
    if file.endswith(end):
        log_files.append(file)
#print(filenames)
for i in range(len(log_files)):
    filename = log_files[i]
    lines_before = []
    #command = "rm " filename.replace("log","gjf")
    #os.system(command)
#open the .com file and read all of it, append it to a list (lines_before)
    with open(filename, 'rb') as f:
        f.seek(-2, 2)
        while f.read(1) != b'\n':
            f.seek(-2, 1)
        last_line = f.readline().decode()
    print(last_line)
    if not "Normal termination" in last_line:
        if version != "":
            command = "obabel -ilog "+ filename + ' -ogjf -O ' + filename[:-4].replace(version,version_new) +'.gjf'
            print(filename)
            os.system(command)
        else:
            command = "obabel -ilog "+ filename + ' -ogjf -O ' + filename[:-4] + version_new +'.gjf'
            print(filename)
            os.system(command)


end2 = version_new + ".gjf"

newlist = []

for file in os.listdir(path):
    if file.endswith(end2):
        newlist.append(file)

print("This is the amount of files to be restarted: " + str(len(newlist)))


print("Now submit those new files\n")

command = "python3 ~/g16_input.py"
os.system(command)

