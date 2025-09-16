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
        lines_before[0]="%chk="+ filename.split("/")[-1][:-4]+"_optint.chk\n"
#change the 2nd line (input line to the one previously defined)
        lines_before[1] = method_line + "\n\n"+ filename[:-4] + "optint" + "\n\n" + " "  + str(charge) + " " + str(multiplicity) + "\n"
        lines_before.append("\n")
        lines_before.append("--Link1--\n")
        lines_before.append("%chk=" + filename.split("/")[-1][:-4] + "_optint.chk\n")
        lines_before.append("# pbe1pbe def2tzvp EmpiricalDispersion=GD3BJ geom=check guess=read scrf(smd, solvent=Dichloromethane) pop=(NBO,CM5)\n")
        lines_before.append("\n")
        lines_before.append("SP\n")
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

def find_ref_coord_indices(filename):
    atoms = ase.io.read(filename)
    #print(atoms)
    indices_to_fix = []
    coords_to_fix = []
    pd_index = None
    # h_index = None
    # c_index = None
    # o_index1 = None
    # o_index2 = None
    for i, atom in enumerate(atoms):
        if atom.symbol == "Pd":
            pd_index = i
            break
    start_index = pd_index - 8
    scaffold_indices = [i for i in range(start_index)]
    #print(scaffold_indices)
    end_index = pd_index + 8
    substrate_indices = [i for i in range(end_index,len(atoms))]

    #print(substrate_indices)
    indices_to_fix = scaffold_indices + substrate_indices 
    for i in substrate_indices:
        coords_to_fix.append(atoms.get_positions()[i])
    return scaffold_indices, substrate_indices, indices_to_fix, coords_to_fix

charge = 0
multiplicity = 1

dir_name = sys.argv[1]
os.chdir(dir_name)
print(dir_name)

#with open('../listofAaronSubs_standard.txt', 'r') as file:
with open('../listofAaronSubs.txt', 'r') as file:
    data = file.readlines()

file_to_sub1 = "Int1.xyz"#sys.argv[2]#"Int1.xyz"
file_to_sub2 = "Int2.xyz"#sys.argv[3]#"Int2.xyz"
atoms_to_subs = sys.argv[2:]#[1,5,8]
print("atoms to substitute:",atoms_to_subs)

Int1_path = "Int1"
Int2_path = "Int2"


try:
    os.system(f"mv {Int2_path} {Int2_path}_old")
except Exception as e:
    print(e)

if not os.path.exists(Int1_path):
    os.mkdir(Int1_path)
if not os.path.exists(Int2_path):
    os.mkdir(Int2_path)

subs_counter = 0
for atom in atoms_to_subs:
    #for i,line in enumerate(data):
    for i in range(len(data)):
        #print(line)
        #print(i)
        line = data[i]
        line = line.strip()[:-4]
        subs_counter += 1
        try:#substitute.py or substitute
            #print("Int1")
            print(subs_counter)
            f_to_sub1 = FileReader(file_to_sub1)
            geom_to_sub1 = Geometry(f_to_sub1)
            #print(file_to_sub)
            #print(geom_to_sub)
            #print("find target:", geom_to_sub1.find(atom))
            #print("sub, target:",line,atom)
            
            sub, sub.atoms, bond, shift, angle, q = geom_to_sub1.substitute(line, atom, minimize=True)
            #print(sub, sub.atoms, bond, shift, angle)
            geom_to_sub1.write(f"{Int1_path}/Int1_Sub{subs_counter}_input")
            create_int_com(f"{Int1_path}/Int1_Sub{subs_counter}_input.xyz")
            os.system(f"mv {Int1_path}/Int1_Sub{subs_counter}_input.com {Int1_path}/Int1_Sub{subs_counter}.com")
            #os.system(f"~/bin/subg16_jed {Int1_path}/Int1_Sub{subs_counter}.com")

            scaffold_indices, substrate_indices, indices_to_fix, coords_to_fix = find_ref_coord_indices(f"{Int1_path}/Int1_Sub{subs_counter}_input.xyz")
            #print(scaffold_indices)
            #print(substrate_indices)
            #print("Int2")
            f_to_sub2 = FileReader(file_to_sub2)
            geom_to_sub2 = Geometry(f_to_sub2)
            
            geom_to_sub2.substitute(line, atom, minimize=False, minimize_RMSD=True,ref_coords=geom_to_sub1.coords, indices_scaffold=scaffold_indices, indices_sub=substrate_indices)#, set_angle=angle)#set_q=q)#, set_shift=shift, set_bond=bond)
            geom_to_sub2.write(f"{Int2_path}/Int2_Sub{subs_counter}_input")
            #print(f"{Int2_path}/Int2_Sub{subs_counter}_input")
            #print(os.getcwd())  
            create_int_com(f"{Int2_path}/Int2_Sub{subs_counter}_input.xyz")
            os.system(f"mv {Int2_path}/Int2_Sub{subs_counter}_input.com {Int2_path}/Int2_Sub{subs_counter}.com")
            #os.system(f"~/bin/subg16_jed {Int2_path}/Int2_Sub{subs_counter}.com")

        except Exception as e:
            print(e) 
