#Juliette Schleicher
#This script automates the Gaussian input files

import glob
import os
import shutil
from pathlib import Path

path = os.getcwd()

charge = 0

multiplicity = 1

submit = str(input("Do you want to submit the jobs right away? y/n \n"))

method = str(input("Which method to use? pbe1pbe as default \n"))

if method == "":
    meth = "pbe1pbe"
else:
    meth = method
 

basis = str(input("Which basis set to use? e.g. def2svp or def2tzvp, default =def2svp \n"))

if basis == "":
    basis = "def2svp"

dispersion = str(input("Dispersion correction? y/n, default=y \n"))

if dispersion == "" or dispersion == "y":
    disp = "EmpiricalDispersion=GD3BJ"
    dis = "D3BJ"
else:
    disp = ""

solvent = str(input("Solvation? y/n \n"))

if solvent == "y":
    solv = "scrf(smd, solvent=Dichloromethane)"
    sol = "SMD"
else:
    solv = ""
    sol = ""

add = input("Anything to add to the filename, e.g. _TZ\n")

filenames_gjf = []

start = input("Starting point: xyz, log or gjf \n")

if start != "gjf":
    newlist = "filenames_" + start
    newlist = []
    for file in os.listdir(path):
        if file.endswith(start):
            newlist.append(file)
#print(filenames)
    for i in range(len(newlist)):
        filename = newlist[i]
        command = 'obabel -i' + start +" "+ filename + ' -ogjf -O ' + filename[:-4] + add + '.gjf'
        print(command)
        os.system(command)



#gets all Gaussian input files, which were created with gview, but without adjusting the method
for file in os.listdir(path):
    if file.endswith(".gjf") and os.path.getsize(file) != 0:
	    filenames_gjf.append(file)
#print(filenames)

#for all the input files (in a dir), the following process is carried out:
for i in range(len(filenames_gjf)):
    filename = filenames_gjf[i]
    
    if "ts" in filename.lower() and "scan" not in filename.lower() and '_f' not in filename.lower() :

        type = "optts"
        method_line = '# ' + meth + " " + basis + " " + disp + " " + solv + " opt=(ts, noeigentest, calcfc, maxstep=10) freq=noraman"
        method_short = meth + "_" + dis + "_" + basis + "_" +  sol + "_optts"
    elif "irc" in filename.lower() and "fw" in filename.lower():
        type = "IRCfw"
        method_line = '# ' + meth + " " + basis + " " + disp + " " + solv + " IRC=(forward,maxpoints=150,recalc=10,calcfc,maxcycle=40,tight,cartesian,lqa,stepsize=15)"
        method_short = meth + "_" + dis + "_" + basis + "_" +  sol + "_IRCfw"
    elif "irc" in filename.lower() and "rv" in filename.lower():
        type = "IRCrv"
        method_line = '# ' + meth + " " + basis + " " + disp + " " + solv + " IRC=(reverse,maxpoints=150,recalc=10,calcfc,maxcycle=40,tight,cartesian,lqa,stepsize=15)"
        method_short = meth + "_" + dis + "_" + basis + "_" +  sol + "_IRCrv"
    elif "int" in filename.lower() and "scan" not in filename.lower() and '_f' not in filename.lower():
        type = "optfreq"
        method_line = '# ' + meth + " " + basis + " " + disp + " " + solv + " opt freq=noraman"
        method_short = meth + "_" + dis + "_" + basis + "_" +  sol + "_optfreq"
    elif "ts" in filename.lower() and "scan" in filename.lower():
        type = "scan"
        method_line = '# ' + meth + " " + basis + " " + disp + " " + solv + " opt=(modredundant)"
        method_short = meth + "_" + dis + "_" + basis + "_" +  sol + "_scan"
    elif "ts" in filename.lower() and "_f" in filename.lower():
        type = "opttsSP"
        method_line = '# ' + meth + " " + basis + " " + disp + " " + solv + " opt=(ts, noeigentest, calcfc, maxstep=10) freq=noraman"
        method_short = meth + "_" + dis + "_" + basis + "_" +  sol + "_opttsSP"
    elif "int" in filename.lower() and "_f" in filename.lower():
        type = "optfreqSP"
        method_line = '# ' + meth + " " + basis + " " + disp + " " + solv + " opt freq=noraman"
        method_short = meth + "_" + dis + "_" + basis + "_" +  sol + "_optfreqSP"
    else:
        type = "sp"
        method_line = '# ' + meth + " " + basis + " " + disp + " " + solv 
        method_short = meth + "_" + dis + "_" + basis + "_" +  sol + "_sp"

    #print(i, filename)
    lines_before = []
#open the .com file and read all of it, append it to a list (lines_before)
    with open(filename, 'r') as f:
        for line in f:
            #print(line)
            lines_before.append(line)
        #print(lines_before)
        
        # Remove connectivity/internal coordinate lines after the geometry block
        # Find the first blank line after the coordinates (line 6+)
        coord_end_idx = -1
        for idx in range(6, len(lines_before)):
            line = lines_before[idx].strip()
            if line == "":
                coord_end_idx = idx
                break
        
        # If we found the blank line after coordinates, remove everything after it
        # except for modredundant statements (lines starting with B, A, D, X)
        if coord_end_idx > 0:
            # Check if there are any modredundant statements after all the connectivity
            modredundant_lines = []
            for idx in range(coord_end_idx + 1, len(lines_before)):
                line = lines_before[idx].strip()
                # Keep only lines that start with B, A, D, or X (modredundant statements)
                if line != "" and len(line) > 0 and line[0] in ['B', 'A', 'D', 'X']:
                    modredundant_lines.append(lines_before[idx])
            
            # Keep everything up to and including the blank line after coords
            # Then add back any modredundant statements
            lines_before = lines_before[:coord_end_idx + 1]
            if modredundant_lines:
                lines_before.extend(modredundant_lines)
                lines_before.append("\n")  # Add final blank line
        
        lines_before[0]="%chk="+ filename[:-4]+".chk\n" 
#change the 2nd line (input line to the one previously defined)
        lines_before[1] = method_line + "\n"
#change the 4th line (title) to the filename without extension
        lines_before[3] = filename[:-4] + method_short + "\n"#"\n"+ filename[:-4] + method_short + "\n"
        lines_before[5] = str(charge) +" " + str(multiplicity) + "\n"
        if type == "IRCfw" or type == "IRCrv":
            lines_before.append("--Link1--\n")
            lines_before.append("%chk=" + filename[:-4] + ".chk\n")
            lines_before.append("# pbe1pbe def2svp EmpiricalDispersion=GD3BJ opt=(maxstep=15) freq=noraman geom=check guess=read\n")
            lines_before.append("\n")
            lines_before.append("Full Opt\n")
            lines_before.append("\n")
            lines_before.append(str(charge) +" " + str(multiplicity) + "\n")
            lines_before.append("\n")
        elif type == "optfreqSP" or type == "opttsSP":
            lines_before.append("--Link1--\n")
            lines_before.append("%chk=" + filename[:-4] + ".chk\n")
            lines_before.append("# pbe1pbe def2tzvp EmpiricalDispersion=GD3BJ scrf(smd, solvent=Dichloromethane) geom=check guess=read pop=(NBO,CM5)\n")
            lines_before.append("\n")
            lines_before.append("SP\n")
            lines_before.append("\n")
            lines_before.append(str(charge) +" " + str(multiplicity) + "\n")
            lines_before.append("\n")
#open the file again, this time with "w" as an option (will overwrite everything, if something goes wrong make sure to have a copy of the files)
    with open(filename[:-4]  + '.com', 'w') as f:
#write the changes
        for element in lines_before:
            f.write(element)
        #print(lines_before)
        
#clean up
 
for i in range(len(filenames_gjf)):
    filename = filenames_gjf[i]
    command = "rm " + filename
    os.system(command)

filenames_com = []

for file in os.listdir(path):
    if file.endswith(".com"):
	    filenames_com.append(file) 
for i in range(len(filenames_com)):
    filename = filenames_com[i]
    isExist = os.path.exists(filename.replace("com","log"))
    #print(isExist)
    if submit == "y":
        if not isExist:
            command = "subg16_jed -p 4 -m 16000 " + filename
            os.system(command)
