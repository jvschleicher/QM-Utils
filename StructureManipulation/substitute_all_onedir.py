import os
import sys
from delete_dup import del_dup 
import ase 
from ase.io import read, write

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


dir_name = sys.argv[1]
os.chdir(dir_name)

with open('../listofAaronSubs.txt', 'r') as file:
    data = file.readlines()


file_to_sub = sys.argv[2]#"Int1.xyz"
atoms_to_subs = sys.argv[3:]#[1,5,8]

#create_xyz_from_gjf(file_to_sub.replace('xyz','gjf'))
#create xyz from .log file
command = f"obabel -ilog ts1_optts_IRCfw.log -oxyz -O ts1_optts_IRCfw.xyz"
os.system(command)
command = f"obabel -ilog ts1_optts_IRCrv.log -oxyz -O ts1_optts_IRCrv.xyz"
os.system(command)

fw = ase.io.read("ts1_optts_IRCfw.xyz")
rv = ase.io.read("ts1_optts_IRCrv.xyz")

# find the palladium index

pd_idx = [i for i, atom in enumerate(fw) if atom.symbol == "Pd"][0] 
print("Pd",pd_idx)
o_idx = pd_idx - 8
print("O",o_idx)
h_idx = pd_idx - 1
print("H",h_idx)

oh_bond_fw = fw.get_distance(o_idx, h_idx)
oh_bond_rv = rv.get_distance(o_idx, h_idx)
#print(oh_bond_fw, oh_bond_rv)

if oh_bond_fw < oh_bond_rv:
    if oh_bond_fw < 1:
        # copy fw to Int2.xyz
        write("Int2.xyz", fw)
        write("Int1.xyz", rv)
    else:
        raise ValueError("The bond length for Int2 is too long, check the structure")
elif oh_bond_rv < oh_bond_fw:
    if oh_bond_rv < 1:
        # copy rv to Int2.xyz
        write("Int2.xyz", rv)
        write("Int1.xyz", fw)
    else:
        raise ValueError("The bond length for Int2 is too long, check the structure")

if not os.path.exists("Int1_v2/"):
    os.mkdir("Int1_v2/")
if not os.path.exists("Int2_v2/"):
    os.mkdir("Int2_v2/")

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
            del_dup(f"{file_to_sub[:-4]}_Sub{subs_counter}.xyz")
            create_int_com(f"{file_to_sub[:-4]}_Sub{subs_counter}.xyz")
            # copy the files to the respective directories
            if "Int1" in file_to_sub:
                os.system(f"cp {file_to_sub[:-4]}_Sub{subs_counter}.com Int1_v2/")
                file_path = f"Int1_v2/{file_to_sub[:-4]}_Sub{subs_counter}.com"
            elif "Int2" in file_to_sub:
                os.system(f"cp {file_to_sub[:-4]}_Sub{subs_counter}.com Int2_v2/")
                file_path = f"Int2_v2/{file_to_sub[:-4]}_Sub{subs_counter}.com"
            command = f"~/bin/subg16 -p 4 -m 16000 {file_path}"
            os.system(command)
        except Exception as e:
            print(e) 
