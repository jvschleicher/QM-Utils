# %%
import os
import numpy as np
# ase
from ase import io
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.geometry.analysis import Analysis as ana
from glob import glob as glob

path = os.getcwd()
os.chdir(path)

# %%
def make_HF3c_input(xyz_file, charge, mult, filename):
    """Generate an input file for HF-3c calculation using ORCA.

    Args:
        xyz file
        charge (int): The total charge of the molecule.
        mult (int): The multiplicity of the molecule.
        filename (str): The name of the output file.
    """

    # read xyz file
    xyz_lines = open(xyz_file).readlines()
    num_atoms = int(xyz_lines[0].strip())
    atoms = []
    coords = []
    for line in xyz_lines[2:2 + num_atoms]:
        parts = line.strip().split()
        atoms.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])


    with open(filename, 'w') as f:
        f.write(f"! HF-3c opt numfreq TightSCF\n")
        f.write(f"* xyz {charge} {mult}\n")
        for atom, coord in zip(atoms, coords):
            f.write(f"  {atom}  {coord[0]:.6f}  {coord[1]:.6f}  {coord[2]:.6f}\n")
        f.write("*\n")

# %%
def make_xtb_Honly(xyz_file, charge, mult, filename):
    """Generate an input file for xTB calculation that only optimised H atoms.

    Args:
        xyz file
        charge (int): The total charge of the molecule.
        mult (int): The multiplicity of the molecule.
        filename (str): The name of the output file.
    """


    # read xyz file
    xyz_lines = open(xyz_file).readlines()
    num_atoms = int(xyz_lines[0].strip())
    atoms = []
    coords = []
    for line in xyz_lines[2:2 + num_atoms]:
        parts = line.strip().split()
        atoms.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])

    with open(filename, 'w') as f:
        f.write(f"! XTB OPT\n")
        f.write(f"%maxcore 1000\n")
        f.write(f"%PAL NPROCS 4 END\n")
        f.write(f"%geom\n")
        f.write(f"optimizehydrogens true\n")
        f.write(f"end\n\n")
        f.write(f"* xyz {charge} {mult}\n")
        for atom, coord in zip(atoms, coords):
            f.write(f"  {atom}  {coord[0]:.6f}  {coord[1]:.6f}  {coord[2]:.6f}\n")
        f.write("*\n")

# %%
def make_xtb_constrained(xyz_file, charge, mult, filename, to_freeze=None):
    """Generate an input file for xTB calculation with constrained atoms .

    Args:
        xyz file
        charge (int): The total charge of the molecule.
        mult (int): The multiplicity of the molecule.
        filename (str): The name of the output file.
    """
   

    # read xyz file
    xyz_lines = open(xyz_file).readlines()
    num_atoms = int(xyz_lines[0].strip())
    atoms = []
    coords = []
    for line in xyz_lines[2:2 + num_atoms]:
        parts = line.strip().split()
        atoms.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])

    if to_freeze is None:
        raise ValueError("needs atoms to freeze!!!")#to_freeze = list(range(0, 86))  # Python 0-based indexing

    with open(filename, 'w') as f:
        f.write(f"! XTB OPT\n")
        f.write(f"%PAL NPROCS 4 END\n")
        f.write(f"%maxcore 1000\n")
        f.write(f"%geom\n")
        f.write(f"  Constraints\n")
        # Each atom needs to be constrained individually with {C atom C}
        for atom_idx in to_freeze:
            f.write(f"    {{C {atom_idx} C}}\n")  # ORCA uses 0-based indexing in constraints
        f.write(f"  end\n")
        f.write(f"end\n\n")
        f.write(f"* xyz {charge} {mult}\n")
        for atom, coord in zip(atoms, coords):
            f.write(f"  {atom}  {coord[0]:.6f}  {coord[1]:.6f}  {coord[2]:.6f}\n")
        f.write("*\n")


# %%
def make_xtb(xyz_file, charge, mult, filename):
    """Generate an input file for xTB calculation with extreme convergence criteria (extreme is not available through orca, use VeryTight instead).

    Args:
        xyz file
        charge (int): The total charge of the molecule.
        mult (int): The multiplicity of the molecule.
        filename (str): The name of the output file.
    """
    

    # read xyz file
    xyz_lines = open(xyz_file).readlines()
    num_atoms = int(xyz_lines[0].strip())
    atoms = []
    coords = []
    for line in xyz_lines[2:2 + num_atoms]:
        parts = line.strip().split()
        atoms.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])

    with open(filename, 'w') as f:
        f.write(f"! XTB VeryTightOPT\n")  
        f.write(f"%PAL NPROCS 4 END\n")
        f.write(f"%maxcore 1000\n")
        f.write(f"%geom\n")
        f.write(f"  MaxIter 500\n")
        #f.write(f"  Convergence VeryTight\n")  # Use VeryTight convergence
        f.write(f"end\n\n")
        f.write(f"* xyz {charge} {mult}\n")
        for atom, coord in zip(atoms, coords):
            f.write(f"  {atom}  {coord[0]:.6f}  {coord[1]:.6f}  {coord[2]:.6f}\n")
        f.write("*\n")
