import os
import sys


def create_int_com(filename):
    add = ""
    method_line = "# pbe1pbe def2svp EmpiricalDispersion=GD3BJ opt freq=noraman"
    lines_before = []
#open the .com file and read all of it, append it to a list (lines_before)
    with open(filename, 'r') as f:
        for line in f:
            #print(line)
            lines_before.append(line)
        #print(lines_before)
        lines_before[0]="%chk="+ filename[:-4]+"_optint.chk\n"
#change the 2nd line (input line to the one previously defined)
        lines_before[1] = method_line + "\n\n"+ filename[:-4] + "optint" + "\n\n" + " "  + str(charge) + " " + str(multiplicity) + "\n"
        lines_before.append("\n")
#open the file again, this time with "w" as an option (will overwrite everything, if something goes wrong make sure to have a copy of the files)
    with open(filename[:-4]  + add + '.com', 'w') as f:
#write the changes
        for element in lines_before:
            f.write(element)

def create_xyz_from_gjf(filename):
    with open(filename,"r") as file:
        lines = file.readlines()
        start = None
        for i, line in enumerate(lines):
            if "0 1" in line:
                start = i+1
                break
            elif "0  1" in line:
                start = i+1
                break
        if start is None:
            raise ValueError("Could not find the coordinates in the gjf file")
    # Extract the coordinate block
    coordinates = []
    for line in lines[start:]:
        if line.strip() == "" or line[0].isdigit():
            break
        coordinates.append(line.strip())    
    #print(coordinates, len(coordinates))    
    with open(filename.replace('gjf','xyz'),"w") as file:
        file.write(f"{len(coordinates)}\n\n")
        for coord in coordinates:
            file.write(coord + "\n")

charge = 0
multiplicity = 1


with open('../listofAaronSubs.txt', 'r') as file:
    data = file.readlines()

file_to_sub = sys.argv[1]#"Int1.xyz"
atoms_to_subs = sys.argv[2:]#[1,5,8]

create_xyz_from_gjf(file_to_sub.replace('xyz','gjf'))

subs_counter = 0
for atom in atoms_to_subs:
    for i,line in enumerate(data):
        #print(line)
        #print(i)
        line = line.strip()[:-4]
        subs_counter += 1
        try:
            command = f"substitute.py {file_to_sub} -s {atom}={line} -m -o {file_to_sub[:-4]}_Sub{subs_counter}.xyz"
            print(command)
            os.system(command)
            create_int_com(f"{file_to_sub[:-4]}_Sub{subs_counter}.xyz")
            command = f"~/bin/subg16 -p 4 -m 16000 {file_to_sub[:-4]}_Sub{subs_counter}.com"
            os.system(command)
        except Exception as e:
            print(e)
