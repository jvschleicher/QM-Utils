# %% [markdown]
# 

# %%
import pandas as pd
import numpy as np
import os
import sys
import ase
from ase.io import read
from glob import glob as glob
# %%

# %%
def fix_atoms(filename):
    # fix all the atoms except the methyl groups from the acetate
    h_index = 20
    ir_index = 19
    c_index = 86
    atoms_to_fix = [h_index, c_index, ir_index]

    return atoms_to_fix
# %%
def make_constrained_ts(filename, add = "_optTS",charge = 0, multiplicity = 1):
    print(filename)
    atoms_to_fix = fix_atoms(filename)
    method_line = "# pbe1pbe def2svp EmpiricalDispersion=GD3BJ opt=(loose,modredundant)"
    lines_before = []
#open the .com file and read all of it, append it to a list (lines_before)
    with open(filename, 'r') as f:
        for line in f:
            #print(line)
            lines_before.append(line)
        #print(lines_before)
        lines_before[0]="%chk="+ filename.split("/")[-1]+"_optts.chk\n"
#change the 2nd line (input line to the one previously defined)
        lines_before[1] = method_line + "\n\n"+ filename[:-4] + "_optts" + "\n\n" + " "  + str(charge) + " " + str(multiplicity) + "\n"
        lines_before.append("\n")
        for atom in atoms_to_fix:
            lines_before.append(f"{atom} F\n")
        lines_before.append("\n")
        lines_before.append("--Link1--\n")
        lines_before.append("%chk=" + filename.split("/")[-1] + "_optts.chk\n")
        lines_before.append("# pbe1pbe def2svp EmpiricalDispersion=GD3BJ geom=check guess=read opt=(ts,noeigentest,calcfc,maxstep=10,nofreeze) freq=noraman\n")
        lines_before.append("\n")
        lines_before.append("Opt\n")
        lines_before.append("\n")
        lines_before.append(str(charge) +" " + str(multiplicity) + "\n")
        lines_before.append("\n")
        lines_before.append("--Link1--\n")
        lines_before.append("%chk=" + filename.split("/")[-1] + "_optts.chk\n")
        lines_before.append("# pbe1pbe def2tzvp EmpiricalDispersion=GD3BJ geom=check guess=read scrf(smd, solvent=Dichloromethane) pop=(NBO,CM5)\n")
        lines_before.append("\n")
        lines_before.append("SP\n")
        lines_before.append("\n")
        lines_before.append(str(charge) +" " + str(multiplicity) + "\n")
        lines_before.append("\n")    
        
    # write to file
    with open(filename.replace('_guessTS.xyz','.com'), 'w') as f:
        for element in lines_before:
            f.write(element)
        #print(lines_before)



# %%
def make_interpolation(dir_name, dir_add):
    # Define the directory paths
    path = os.getcwd()
    #print(os.getcwd())
    if dir_name != "":
        os.chdir(f"{path}/{dir_name}")

        int1_dir = os.path.join(f"{dir_name}", f'Int1{dir_add}')
        int2_dir = os.path.join(f"{dir_name}", f'Int2{dir_add}')
        qst2_dirname = f"interpolation{dir_add}2"
        # Create a dictionary to make the qst2 input files in 
        if not os.path.exists(f"{path}/{dir_name}/{qst2_dirname}"):
            os.mkdir(f"{path}/{dir_name}/{qst2_dirname}")
    else:
        os.chdir(f"{path}")

        int1_dir = os.path.join(f"{path}", f'Int1{dir_add}')
        int2_dir = os.path.join(f"{path}", f'Int2{dir_add}')
        qst2_dirname = f"interpolation{dir_add}2"
        # Create a dictionary to make the qst2 input files in
        if not os.path.exists(f"{path}/{qst2_dirname}"):
            os.mkdir(f"{path}/{qst2_dirname}")


    #os.chdir(int1_dir)
    #os.system("bash /home/jschleic/scripts/xyz/makexyz.sh")
    #os.chdir(int2_dir)
    #os.system("bash /home/jschleic/scripts/xyz/makexyz.sh")
    #os.chdir(f"{base}/{dir_name}")
    # Get the list of files in the directories
    #int1_files = os.listdir(f"{path}/{int1_dir}")
    #int2_files = os.listdir(f"{path}/{int2_dir}")
    int1_files = os.listdir(f"{int1_dir}")
    int2_files = os.listdir(f"{int2_dir}")

    #print(int1_files)
    # Filter the files based on the conditions
    intermediates = []
    for file in int1_files:
        if file.endswith('.log'):
            #print(file)
            #print(file.split("1_")[-1])
            sub_no1 = file.split("Int1_")[-1]
            for file in int2_files:
                if file.endswith('.log'):
                    #print(file)
                    #print(file.split("2_")[-1])
                    sub_no2 = file.split("Int2_")[-1]
                    if sub_no1 == sub_no2:
                        intermediates.append(sub_no1)


    # Sort the intermediates
    intermediates.sort()


    # %%
    # Print the ordered intermediates
    #for i, intermediate in enumerate(intermediates):
    #    print(f"sub{i+1}: {intermediate}")


    # %%
    # Create a dictionary to make the qst2 input files in 

    for i, intermediate in enumerate(intermediates):
    #for i in range(10):
        #intermediate = intermediates[i]
        print(intermediate)
        #print(f"{int1_dir.split('/')[-1]}/Int1_*{intermediate}")
        int1_log_name = glob(f"{int1_dir.split('/')[-1]}/Int1_*{intermediate}")[0]
    
        #print(f"{int1_dir.split('/')[-1]}/Int1_*{intermediate}")

        int2_log_name = glob(f"{int2_dir.split('/')[-1]}/Int2_*{intermediate}")[0]
        int1_xyz = read(int1_log_name)    #f"{int1_dir.split('/')[-1]}/Int1_{intermediate}"
        int2_xyz = read(int2_log_name)   
        # get the atoms and positions of the intermediates, same format as in the xyz files
        atoms1 = int1_xyz.get_chemical_symbols()
        positions1 = int1_xyz.get_positions()
        atoms2 = int2_xyz.get_chemical_symbols()
        positions2 = int2_xyz.get_positions()
        #get scf energy of the intermediates
        try:
            scf1 = int1_xyz.get_potential_energy()
            scf2 = int2_xyz.get_potential_energy()
        except Exception as e:
            print(f"Error: {e}")
            scf1 = None
            scf2 = None
        # create the qst2 input file
        #with open(f"{path}/{dir_name}/{qst2_dirname}/interpolation_{intermediate.replace('.log','_input.xyz')}", 'w') as f:
        with open(f"{path}/{qst2_dirname}/interpolation_{intermediate.replace('.log','_input.xyz')}", 'w') as f:

            f.write(f"{len(atoms1)}\n")
            f.write(f"{scf1}\n")
            for atom, position in zip(atoms1, positions1):
                f.write(f"{atom} {position[0]} {position[1]} {position[2]}\n")
            f.write(f"{len(atoms2)}\n")
            f.write(f"{scf2}\n")

            for atom, position in zip(atoms2, positions2):
                f.write(f"{atom} {position[0]} {position[1]} {position[2]}\n")
        #print(f"python -m geodesic_interpolate {path}/{dir_name}/{qst2_dirname}/interpolation_{intermediate.replace('log','xyz')} --nimages 3 --output {path}/{dir_name}/{qst2_dirname}/interpolation_{intermediate.replace('.log','_out.xyz')}")
        os.system(f"python -m geodesic_interpolate {path}/{qst2_dirname}/interpolation_{intermediate.replace('.log','_input.xyz')} --nimages 3 --output {path}/{qst2_dirname}/interpolation_{intermediate.replace('.log','_out.xyz')}")

        # read the output file, return the second xyz structure
        #out = ase.io.read(f"{path}/{qst2_dirname}/interpolation_{intermediate.replace('.log','_out.xyz')}", index=1)
        #atoms_out = out.get_chemical_symbols()
        #positions_out = out.get_positions()
        #with open(f"{path}/{qst2_dirname}/interpolation_{intermediate.replace('.log','_guessTS.xyz')}", 'w') as f:
        #    f.write(f"{len(atoms_out)}\n")
        #    f.write(f"\n")
        #    for atom, position in zip(atoms_out, positions_out):
        #        f.write(f"{atom} {position[0]} {position[1]} {position[2]}\n")
        #make_constrained_ts(f"{path}/{qst2_dirname}/interpolation_{intermediate.replace('.log','_guessTS.xyz')}")
        #submit
        #os.system(f"~/bin/subg16 {path}/{dir_name}/{qst2_dirname}/interpolation_{intermediate.replace('.log','.com')}") 
        #os.system(f"mv {path}/{dir_name}/{qst2_dirname}/interpolation_{intermediate.replace('.log','_guessTS.xyz')[:-4]}.com {path}/{dir_name}/{qst2_dirname}/interpolation_{intermediate.replace('.log','.com')}")
    os.chdir(path)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        make_interpolation(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 2:
        make_interpolation(sys.argv[1], "")
    else:
        make_interpolation("","")
