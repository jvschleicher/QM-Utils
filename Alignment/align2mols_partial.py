
#!/usr/bin/env python3
"""
Align two molecular structure files using the Kabsch algorithm,
optionally on a subset of atoms.

Reads any format supported by ASE (XYZ, Gaussian log, ORCA output, etc.).
The second file is aligned onto the first (reference).

Outputs:
  <file1_stem>_aligned_with_<file2_stem>.xyz  – file1 in its original position
  <file2_stem>_aligned_with_<file1_stem>.xyz  – file2 rotated/translated onto file1
  aligned_<file1_stem>_<file2_stem>.xyz        – both structures in one XYZ file

Usage:
  python align2mols_partial.py file1.xyz file2.xyz --atoms 1 2 5 6
  python align2mols_partial.py file1.xyz file2.xyz          # use all atoms
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


def align_to_reference(atoms_ref, atoms_to_align, align_indices=None):
    """
    Align atoms_to_align to atoms_ref using specified atom indices.
    If align_indices is None, use all atoms for alignment.
    
    Args:
        atoms_ref: ASE Atoms object (reference structure)
        atoms_to_align: ASE Atoms object (structure to align)
        align_indices: List of atom indices (0-indexed) to use for alignment.
                       If None, all atoms are used.
    
    Returns:
        ASE Atoms object with aligned coordinates
    """
    # Get coordinates
    P = atoms_ref.get_positions()
    Q = atoms_to_align.get_positions()
    
    # If no alignment indices specified, use all atoms
    if align_indices is None:
        align_indices = list(range(len(P)))
    
    # Extract only the atoms to align on
    P_align = P[align_indices]
    Q_align = Q[align_indices]
    
    # Center the alignment atoms
    P_centroid = centroid(P_align)
    Q_centroid = centroid(Q_align)
    
    P_align_centered = P_align - P_centroid
    Q_align_centered = Q_align - Q_centroid
    
    # Calculate rotation using only the alignment atoms
    U = kabsch(P_align_centered, Q_align_centered)
    
    # Apply rotation to all atoms
    # First center all atoms using Q's centroid
    Q_centered = Q - Q_centroid
    # Rotate all atoms
    Q_rotated = np.dot(Q_centered, U.T)
    # Translate to P's centroid
    Q_final = Q_rotated + P_centroid
    
    # Create new Atoms object with aligned coordinates
    atoms_aligned = atoms_to_align.copy()
    atoms_aligned.set_positions(Q_final)

    return atoms_aligned


def main():
    parser = argparse.ArgumentParser(
        description='Align two molecular structures using the Kabsch algorithm (full or partial).'
    )
    parser.add_argument('file1', help='Reference structure (any ASE-readable format)')
    parser.add_argument('file2', help='Structure to align onto file1 (any ASE-readable format)')
    parser.add_argument(
        '--atoms', '-a',
        nargs='+',
        type=int,
        default=None,
        metavar='N',
        help='1-based atom indices to use for alignment (e.g. --atoms 1 2 5 6). '
             'If omitted, all atoms are used.'
    )
    parser.add_argument(
        '--index', '-i',
        default='-1',
        help='Frame index to read from each file (default: -1, i.e. last frame). '
             'Use "0" for first frame.'
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

    # Convert 1-based CLI indices to 0-based
    align_indices = None
    if args.atoms is not None:
        align_indices = [i - 1 for i in args.atoms]
        out_of_range = [i for i in align_indices if not (0 <= i < len(atoms_ref))]
        if out_of_range:
            raise ValueError(
                f"Atom indices out of range (structure has {len(atoms_ref)} atoms): "
                f"{[i + 1 for i in out_of_range]}"
            )
        print(f"Aligning on {len(align_indices)} atoms: {args.atoms}")
    else:
        print(f"Aligning on all {len(atoms_ref)} atoms")

    atoms_aligned = align_to_reference(atoms_ref, atoms_mob, align_indices)

    stem1 = path1.stem
    stem2 = path2.stem

    out_ref      = path1.parent / f"{stem1}_aligned_with_{stem2}_on_{'_'.join(map(str, args.atoms))}.xyz" if args.atoms else path1.parent / f"{stem1}_aligned_with_{stem2}.xyz"
    out_mob      = path2.parent / f"{stem2}_aligned_with_{stem1}_on_{'_'.join(map(str, args.atoms))}.xyz" if args.atoms else path2.parent / f"{stem2}_aligned_with_{stem1}.xyz"
    out_combined = path1.parent / f"aligned_{stem1}_{stem2}.xyz"

    write(str(out_ref),      atoms_ref,     format='xyz')
    write(str(out_mob),      atoms_aligned, format='xyz')
    write(str(out_combined), [atoms_ref, atoms_aligned], format='xyz')

    print(f"Reference (unchanged): {out_ref}")
    print(f"Aligned structure:     {out_mob}")
    print(f"Combined file:         {out_combined}")


if __name__ == "__main__":
    main()