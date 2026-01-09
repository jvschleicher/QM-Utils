from ase.io import read, write
import sys

file = sys.argv[1]

import cclib

#data = cclib.io.ccread(file)
#print(len(data.atomcoords))  # number of optimization steps
from cclib.io import xyzwriter

#xyzwriter.XYZ(data).write(f"{file.replace('.out','_traj.xyz')}")


from ase.io import read, write

#atoms = read(file, index=":")


import sys

outfile = sys.argv[1]
frames = []
current = []

with open(outfile) as f:
    lines = iter(f)
    for line in lines:
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
            current = []
            next(lines)  # dashed line
            # Don't use continue - process atoms immediately after this
            for atom_line in lines:
                if atom_line.strip() == "":
                    if current:
                        frames.append(current)
                    current = None
                    break
                
                parts = atom_line.split()
                if len(parts) == 4:
                    current.append(parts)

if not frames:
    print("No ANGSTROEM geometries found.")
    sys.exit(1)

with open(f"{file.replace('.out','_traj.xyz')}", "w") as out:
    for i, frame in enumerate(frames):
        out.write(f"{len(frame)}\n")
        out.write(f"ORCA optimization step {i}\n")
        for el, x, y, z in frame:
            out.write(f"{el:2s} {x} {y} {z}\n")




