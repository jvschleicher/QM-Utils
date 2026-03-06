#!/usr/bin/env python3
"""
Align two molecular structure files using the Kabsch algorithm. 
Fully based on https://github.com/charnley/rmsd/tree/master 

Reads any format supported by ASE (XYZ, Gaussian log, ORCA output, etc.).
The second file is aligned onto the first (reference).

Outputs:
  <file1_stem>_aligned_with_<file2_stem>.xyz  – file1 in its original position
  <file2_stem>_aligned_with_<file1_stem>.xyz  – file2 rotated/translated onto file1
  aligned_<file1_stem>_<file2_stem>.xyz        – both structures in one XYZ file

Usage:
  python align2mols.py file1.xyz file2.xyz
"""

import argparse
import numpy as np
from ase.io import read, write
from pathlib import Path


def centroid(X):
    """Calculate the centroid of a set of coordinates."""
    return X.mean(axis=0)


def kabsch(P, Q):
    """
    Calculate the optimal rotation matrix using Kabsch algorithm.
    
    Args:
        P: Nx3 matrix of coordinates (reference)
        Q: Nx3 matrix of coordinates (to be aligned)
    
    Returns:
        U: 3x3 rotation matrix
    """
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    return U


def align_2_mol(atoms_ref, atoms_to_align):
    """
    Align two molecules using all atoms.
    
    Args:
        atoms_ref: ASE Atoms object (reference structure)
        atoms_to_align: ASE Atoms object (structure to align)
    
    Returns:
        ASE Atoms object with aligned coordinates
    """
    P = atoms_ref.get_positions()
    Q = atoms_to_align.get_positions()
    
    P_centroid = centroid(P)
    Q_centroid = centroid(Q)
    
    P_centered = P - P_centroid
    Q_centered = Q - Q_centroid
    
    U = kabsch(P_centered, Q_centered)
    
    Q_rotated = np.dot(Q_centered, U.T)
    Q_final = Q_rotated + P_centroid
    
    atoms_aligned = atoms_to_align.copy()
    atoms_aligned.set_positions(Q_final)
    
    return atoms_aligned


def main():
    parser = argparse.ArgumentParser(
        description='Align two molecular structures using the Kabsch algorithm.'
    )
    parser.add_argument('file1', help='Reference structure (any ASE-readable format)')
    parser.add_argument('file2', help='Structure to align onto file1 (any ASE-readable format)')
    parser.add_argument(
        '--index', '-i',
        default='-1',
        help='Frame index to read from each file (default: -1, i.e. last frame). '
             'Use "0" for first frame, ":" for all frames.'
    )
    args = parser.parse_args()

    path1 = Path(args.file1)
    path2 = Path(args.file2)

    atoms_ref = read(str(path1), index=args.index)
    atoms_mob = read(str(path2), index=args.index)

    # read() with a scalar index returns a single Atoms object;
    # with ':' it returns a list — normalise to single object
    if isinstance(atoms_ref, list):
        atoms_ref = atoms_ref[-1]
    if isinstance(atoms_mob, list):
        atoms_mob = atoms_mob[-1]

    if len(atoms_ref) != len(atoms_mob):
        raise ValueError(
            f"Atom count mismatch: {path1.name} has {len(atoms_ref)} atoms, "
            f"{path2.name} has {len(atoms_mob)} atoms."
        )

    atoms_aligned = align_2_mol(atoms_ref, atoms_mob)

    stem1 = path1.stem
    stem2 = path2.stem

    out_ref     = path1.parent / f"{stem1}_aligned_with_{stem2}.xyz"
    out_mob     = path2.parent / f"{stem2}_aligned_with_{stem1}.xyz"
    out_combined = path1.parent / f"aligned_{stem1}_{stem2}.xyz"

    write(str(out_ref),      atoms_ref,     format='xyz')
    write(str(out_mob),      atoms_aligned, format='xyz')
    write(str(out_combined), [atoms_ref, atoms_aligned], format='xyz')

    print(f"Reference (unchanged): {out_ref}")
    print(f"Aligned structure:     {out_mob}")
    print(f"Combined file:         {out_combined}")


if __name__ == "__main__":
    main()
