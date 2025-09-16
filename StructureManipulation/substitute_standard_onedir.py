import os
import sys
from delete_dup import del_dup
import ase 
from ase.io import read, write
from AaronTools.fileIO import FileReader, read_types
from AaronTools.geometry import Geometry
from AaronTools.substituent import Substituent


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
        lines_before.append("--Link1--\n")
        lines_before.append("%chk=" + filename[:-4] + "_optint.chk\n")
        lines_before.append("# pbe1pbe def2tzvp EmpiricalDispersion=GD3BJ geom=check guess=read scrf(smd, solvent=Dichloromethane) pop=(NBO,CM5)\n")
        lines_before.append("\n")
        lines_before.append("Full Opt\n")
        lines_before.append("\n")
        lines_before.append(str(charge) +" " + str(multiplicity) + "\n")
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

dir_name = sys.argv[1]
os.chdir(dir_name)

with open('../listofAaronSubs_standard.txt', 'r') as file:
#with open('../listofAaronSubs.txt', 'r') as file:
    data = file.readlines()

file_to_sub1 = sys.argv[2]#"Int1.xyz"
file_to_sub2 = sys.argv[3]#"Int2.xyz"
atoms_to_subs = sys.argv[4:]#[1,5,8]


Int1_path = "Int1"
Int2_path = "Int2"

if not os.path.exists(Int1_path):
    os.mkdir(Int1_path)
if not os.path.exists(Int2_path):
    os.mkdir(Int2_path)

subs_counter = 0
for atom in atoms_to_subs:
    #for i,line in enumerate(data):
    for i in range(len(data))[:5]:
        #print(line)
        #print(i)
        line = data[i]
        line = line.strip()[:-4]
        subs_counter += 1
        try:#substitute.py or substitute
            print("Int1")
            f_to_sub1 = FileReader(file_to_sub1)
            geom_to_sub1 = Geometry(f_to_sub1)
            #print(file_to_sub)
            #print(geom_to_sub)
            #print("find target:", geom_to_sub1.find(atom))
            print("sub, target:",line,atom)
            geom_to_sub1.substitute(line, int(atom), minimize=True)
            print("Int2")
            f_to_sub2 = FileReader(file_to_sub2)
            geom_to_sub2 = Geometry(f_to_sub2)

            geom_to_sub2.substitute(line, int(atom), minimize=True)
            #command = f"substitute {file_to_sub} -s {atom}={line} -m -o {file_to_sub[:-4]}_Sub{subs_counter}.xyz"#-m for minimise -m -o
            #print(command)
            #os.system(command)
            #del_dup(f"{file_to_sub[:-4]}_Sub{subs_counter}.xyz")
            #create_int_com(f"{file_to_sub[:-4]}_Sub{subs_counter}.xyz")
            # copy the files to the respective directories
            if "Int1" in file_to_sub:
                #os.system(f"mv {file_to_sub[:-4]}_Sub{subs_counter}.* {Int1_path}/")
                #file_path = f"{Int1_path}/{file_to_sub[:-4]}_Sub{subs_counter}.com"
                #command = f"~/bin/subg16_jed -p 4 -m 16000 {file_path}"
                #os.system(command)
                pass
            elif "Int2" in file_to_sub:
                os.system(f"mv {file_to_sub[:-4]}_Sub{subs_counter}.* {Int2_path}/")
                file_path = f"{Int2_path}/{file_to_sub[:-4]}_Sub{subs_counter}.com"
                #command = f"~/bin/subg16_jed -p 4 -m 16000 {file_path}"
                #os.system(command)
        except Exception as e:
            print(e) 
