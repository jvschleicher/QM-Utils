import kallisto
from glob import glob as glob
import os
import shutil
import sys

from  AaronTools.geometry import Geometry
from AaronTools.substituent import Substituent

# ase
from ase.io import read
from ase.io import write
from ase.geometry.analysis import Analysis as ana
from ase.neighborlist import *
from ase import Atoms

import re
import pandas as pd

def read_xyz(filename):
    with open(filename, 'r') as file:
        num_atoms = int(file.readline())
        xyz_data = [file.readline() for _ in range(num_atoms+1)]
    return xyz_data

def write_xyz_file(filename, xyz_data):
    with open(filename, 'w') as file:
        for i,line in enumerate(xyz_data):
            file.write(f"{line}")

def extract_xyz_fragments_benz(xyz_data):
    fragment1 = []
    fragment2 = []
    fragment3 = []
    fragment4 = []
    # Find the index where to split the XYZ data
    pd_index = 0
    for i, line in enumerate(xyz_data):
        if line.startswith('Pd'):
            pd_index = i
            break
    # Extract fragment 1
    fragment1 = xyz_data[pd_index+8:]
    
    # Extract fragment 2
    h_index = 0
    for i, line in enumerate(xyz_data):
        if line.startswith('H'):
            h_index = i
            break
    fragment2 = [xyz_data[h_index]]

    # Extract fragment 3
    fragment3 = xyz_data[h_index+1:pd_index]

    # Extract fragment 4
    fragment4 = xyz_data[pd_index:pd_index+8]
    pd_idx = len(fragment1)+2
    info = [f"{len(fragment1)+len(fragment2)+len(fragment3)+len(fragment4)}\n",f"F:1-{pd_idx-1+7};1-{pd_idx+7};{pd_idx-1+7}-{pd_idx+7};\n"]

    return info,fragment1, fragment2, fragment3, fragment4

def extract_xyz_fragments(xyz_data):
    fragment1 = []
    fragment2 = []
    fragment3 = []
    fragment4 = []
    # Find the index where to split the XYZ data
    pd_index = 0
    for i, line in enumerate(xyz_data):
        if line.startswith('Pd'):
            pd_index = i
            break
    # Extract fragment 1 (substrate)
    fragment1 = xyz_data[pd_index+8:]
    
    # Extract fragment 2 (hydrogen)
    h_index = pd_index -1 
    #for i, line in enumerate(xyz_data):
    #    if line.startswith('H'):
    #        h_index = i
    #        
    #        break
    print("h",h_index)
    fragment2 = [xyz_data[h_index]]

    # Extract fragment 3 (second acetate, half bound)
    fragment3 = xyz_data[1:h_index]#xyz_data[h_index+1:pd_index]#get rid of the created with callisto line

    # Extract fragment 4 (first acetate, fully bound)
    fragment4 = xyz_data[pd_index:pd_index+8]
    pd_idx = len(fragment1)+2
    info = [f"{len(fragment1)+len(fragment2)+len(fragment3)+len(fragment4)}\n",f"F:1-{pd_idx-1+7};1-{pd_idx+7};{pd_idx-1+7}-{pd_idx+7};\n"]

    return info,fragment1, fragment2, fragment3, fragment4

#mol_xyz_sp2 = read("ts1_benz.xyz")
#pd_indices_sp2 = [i for i, atom in enumerate(mol_xyz_sp2) if atom.symbol == 'Pd']
#pd_centre_sp2 = pd_indices_sp2[0]
mol_xyz_sp3 = read("ts1_eth.xyz")
pd_indices_sp3 = [i for i, atom in enumerate(mol_xyz_sp3) if atom.symbol == 'Pd']
pd_centre_sp3 = pd_indices_sp3[0]

# get list of dirs that do not have _H in them
dirs = [sys.argv[1]]#[d for d in os.listdir() if os.path.isdir(d) and "_H" not in d]
base_path = os.getcwd()
# then do this for all of the dirs
for d in dirs:

    os.chdir(base_path)
    #shutil.copy("ts1_benz.xyz",f"{d}/ts1_benz.xyz")
    shutil.copy("ts1_eth.xyz",f"{d}/ts1_eth.xyz")
    print(d)
    os.chdir(d)

    path = os.getcwd()

    #turn gjf file into xyz - get the coordinates found beneath charge and multiplicity ( 0 1)
    # get the coordinates and write them to a new file
    with open("subs.gjf","r") as file:
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
    with open("subs.xyz","w") as file:
        file.write(f"{len(coordinates)}\n\n")
        for coord in coordinates:
            file.write(coord + "\n")


    filenames = glob("sub*xyz")
    for i in range(len(filenames)):
        
        Cindices = []
        Hindices = []
        if "plusPh" not in filenames[i]:
            os.chdir(path) # make sure to return to the submit dir
            filename = filenames[i]
            conf_number = filename.split("conformer")[-1].strip(".xyz")
            if not os.path.exists(f"confs{conf_number}"):
                os.mkdir(f"confs{conf_number}")
            
            shutil.copy(filename, f"confs{conf_number}/{filename}")
            #shutil.copy("ts1_benz.xyz",f"confs{conf_number}/ts1_benz.xyz")
            shutil.copy("ts1_eth.xyz",f"confs{conf_number}/ts1_eth.xyz")
            os.chdir(f"confs{conf_number}") # change dir so that all of this script happens for one conformer at the time without overwriting anything
    
            mol_xyz = read(filename)
        #
            analysis = ana(mol_xyz)
            CHBonds = analysis.get_bonds("C", "H", unique=True)
            # CHBonds is a list of a list of tuples
            print(filename, "ch bonds", CHBonds[0])
            for j in range(len(CHBonds[0])):
                Cidx = CHBonds[0][j][0]
                Hidx = CHBonds[0][j][1]
                Cindices.append(Cidx)
                Hindices.append(Hidx)


            for i, (Cidx,Hidx) in enumerate(zip(Cindices,Hindices)):
                mol_xyz = read(filename) # redo the read part to re-initialise it
                # Define a cutoff distance for neighbor detection
                cutoff_distance = 1.0  # Adjust as neede

                # Create a NeighborList object
                nl = NeighborList([cutoff_distance / 2.0] * len(mol_xyz), self_interaction=False, bothways=True)

                # Update the NeighborList with the positions of atoms in the molecule
                nl.update(mol_xyz)
                indices = []
                indices, offsets = nl.get_neighbors(Cidx)
                num_neighbors = len(indices)
                print(filename, Cidx, indices)
                #print("Neighbors",num_neighbors)
                mol_xyz.pop(Hidx)
                filename_pop = f"{filename[:-4]}_pop{Hidx}.xyz"
                write(filename_pop, mol_xyz)
                try:
                    command = f"kallisto sort --start {Cidx} --out {filename[:-4]}_{Cidx + 1}_{Hidx +1}_sorted.xyz {filename_pop}"
                    #print(command)
                    os.system(command)
                except Exception as e:
                    print(f'Exception {e} for {filename}')

                print("created the sorted substrate-structure")


                try:
                    #if num_neighbors <= 3:
                    #    command1 = f"kallisto exs --center {pd_centre_sp2} --subnr 0 ts1_benz.xyz {filename[:-4]}_{Cidx + 1}_{Hidx +1}_sorted.xyz"
                    #    print(command1)
                    #    os.system(command1)
                    #elif num_neighbors == 4:
                    command1 = f"kallisto exs --center {pd_centre_sp3} --subnr 0 ts1_eth.xyz {filename[:-4]}_{Cidx + 1}_{Hidx +1}_sorted.xyz"
                    #print(command1)
                    os.system(command1)

                    command_opt = "xtb newstructure.xyz --opt tight --input constrain.inp"
                    #print(command_opt)
                    #os.system(command_opt)
                    shutil.move("newstructure.xyz",f"ts_{Cidx + 1}_{Hidx +1}.xyz")
                except Exception as e:
                    print(f'Exception {e} for {filename}')

                target_dir = f"{path}_H{Hidx +1}/"
                if not os.path.exists(target_dir):
                    print(target_dir)
                    os.mkdir(target_dir)
                xyz_data = read_xyz(f"ts_{Cidx + 1}_{Hidx +1}.xyz")
                info, fragment1, fragment2, fragment3, fragment4 = extract_xyz_fragments(xyz_data)
                xyz_data = info + fragment1 + fragment3 + fragment2 + fragment4
                #print(info, fragment1, fragment2, fragment3, fragment4)
                write_xyz_file(f"{target_dir}ts1.xyz",xyz_data)
 
 
        #    with open(f"ts_{Cidx + 1}_{Hidx +1}.xyz","r") as file:
        #        for line in file:
        #            print(line)
        #            lines_before.append(line)
        ####continue here
