#!/usr/bin/env python3
"""
Convert Gaussian .com files to .xyz format
"""

import sys
import os
import re


def parse_com_file(com_file):
    """
    Parse a Gaussian .com file and extract molecular coordinates.
    
    Args:
        com_file (str): Path to the .com file
        
    Returns:
        tuple: (atoms, coords, title) where atoms is a list of element symbols,
               coords is a list of (x, y, z) tuples, and title is the molecule title
    """
    atoms = []
    coords = []
    title = ""
    
    with open(com_file, 'r') as f:
        lines = f.readlines()
    
    # Find the geometry section
    # Skip lines until we find charge and multiplicity
    geometry_start = -1
    for i, line in enumerate(lines):
        # Look for title section (after route section)
        if i > 0 and line.strip() and not line.startswith('%') and not line.startswith('#'):
            if not title and i > 1:
                title = line.strip()
        
        # Look for charge/multiplicity line (e.g., "0 1")
        if re.match(r'^\s*[+-]?\d+\s+\d+\s*$', line):
            geometry_start = i + 1
            break
    
    if geometry_start == -1:
        raise ValueError("Could not find geometry section in .com file")
    
    # Parse geometry lines
    for i in range(geometry_start, len(lines)):
        line = lines[i].strip()
        
        # Stop at blank line or end of geometry
        if not line:
            break
            
        # Parse coordinate line: Element X Y Z
        parts = line.split()
        if len(parts) >= 4:
            element = parts[0]
            try:
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                atoms.append(element)
                coords.append((x, y, z))
            except (ValueError, IndexError):
                # Not a valid coordinate line, might be end of geometry
                break
    
    if not atoms:
        raise ValueError("No atoms found in .com file")
    
    return atoms, coords, title


def write_xyz_file(xyz_file, atoms, coords, title=""):
    """
    Write molecular coordinates to .xyz file.
    
    Args:
        xyz_file (str): Path to output .xyz file
        atoms (list): List of element symbols
        coords (list): List of (x, y, z) coordinate tuples
        title (str): Title/comment line
    """
    with open(xyz_file, 'w') as f:
        # Write number of atoms
        f.write(f"{len(atoms)}\n")
        
        # Write title/comment line
        f.write(f"{title if title else 'Converted from .com file'}\n")
        
        # Write coordinates
        for atom, (x, y, z) in zip(atoms, coords):
            f.write(f"{atom:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")


def convert_com_to_xyz(com_file, xyz_file=None):
    """
    Convert a .com file to .xyz format.
    
    Args:
        com_file (str): Path to input .com file
        xyz_file (str): Path to output .xyz file (optional, default: same name with .xyz extension)
    """
    # Generate output filename if not provided
    if xyz_file is None:
        xyz_file = os.path.splitext(com_file)[0] + '.xyz'
    
    # Parse .com file
    atoms, coords, title = parse_com_file(com_file)
    
    # Write .xyz file
    write_xyz_file(xyz_file, atoms, coords, title)
    
    print(f"Converted {com_file} to {xyz_file}")
    print(f"  Atoms: {len(atoms)}")


def main():
    """Main function to handle command line arguments."""
    if len(sys.argv) < 2:
        print("Usage: python com_to_xyz.py <input.com> [output.xyz]")
        print("  If output.xyz is not specified, it will use the same name as input with .xyz extension")
        sys.exit(1)
    
    com_file = sys.argv[1]
    xyz_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.exists(com_file):
        print(f"Error: File '{com_file}' not found")
        sys.exit(1)
    
    try:
        convert_com_to_xyz(com_file, xyz_file)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
