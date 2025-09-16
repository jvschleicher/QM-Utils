#!/usr/bin/env python

import re
import sys
from os.path import dirname, splitext

import networkx as nx
import numpy as np
import scipy.spatial
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


ha_to_kcalmol = 627.509474
kcalmol_to_ha = 0.00159360
symbol_to_number = {
    "Em": 0,  # empty site
    "Vc": 0,  # empty site
    "Va": 0,  # empty site
    "D": 1,
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Me": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Uun": 110,
    "Uuu": 111,
    "Uub": 112,
}

number_to_symbol = {v: k for k, v in symbol_to_number.items()}

symbol_to_mass = {
    "Em": 0.000,
    "Va": 0.000,
    "Vc": 0.000,
    "H": 1.0075,
    "D": 2.01410178,
    "He": 4.002,
    "Li": 6.9675,
    "Be": 9.012,
    "B": 10.8135,
    "C": 12.0106,
    "N": 14.0065,
    "O": 15.999,
    "F": 18.998,
    "Ne": 20.1797,
    "Na": 22.989,
    "Mg": 24.3050,
    "Al": 26.981,
    "Si": 28.085,
    "P": 30.973,
    "S": 32.0675,
    "Cl": 35.4515,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938,
    "Fe": 55.845,
    "Co": 58.933,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.63,
    "As": 74.921,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.905,
    "Zr": 91.224,
    "Nb": 92.906,
    "Mo": 95.96,
    "Tc": 98,
    "Ru": 101.07,
    "Rh": 102.905,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.60,
    "I": 126.904,
    "Xe": 131.293,
    "Cs": 132.905,
    "Ba": 137.327,
    "La": 138.905,
    "Ce": 140.116,
    "Pr": 140.907,
    "Nd": 144.242,
    "Pm": 145,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.925,
    "Dy": 162.500,
    "Ho": 164.930,
    "Er": 167.259,
    "Tm": 168.934,
    "Yb": 173.054,
    "Lu": 174.9668,
    "Hf": 178.49,
    "Ta": 180.947,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.966,
    "Hg": 200.59,
    "Tl": 204.3835,
    "Pb": 207.2,
    "Bi": 208.980,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226,
    "Ac": 227,
    "Th": 232.038,
    "Pa": 231.035,
    "U": 238.028,
    "Np": 237,
    "Pu": 244,
    "Am": 243,
    "Cm": 247,
    "Bk": 247,
    "Cf": 251,
    "Es": 252,
    "Fm": 257,
    "Md": 258,
    "No": 259,
    "Lr": 262,
    "Rf": 267,
    "Db": 268,
    "Sg": 271,
    "Bh": 272,
    "Hs": 270,
    "Mt": 276,
    "Ds": 281,
    "Rg": 280,
    "Cn": 285,
}


missing = 0.2
covalent_radii = np.array(
    [
        missing,  # X
        0.31,  # H
        0.28,  # He
        1.28,  # Li
        0.96,  # Be
        0.84,  # B
        0.76,  # C
        0.71,  # N
        0.66,  # O
        0.57,  # F
        0.58,  # Ne
        1.66,  # Na
        1.41,  # Mg
        1.21,  # Al
        1.11,  # Si
        1.07,  # P
        1.05,  # S
        1.02,  # Cl
        1.06,  # Ar
        2.03,  # K
        1.76,  # Ca
        1.70,  # Sc
        1.60,  # Ti
        1.53,  # V
        1.39,  # Cr
        1.39,  # Mn
        1.32,  # Fe
        1.26,  # Co
        1.24,  # Ni
        1.32,  # Cu
        1.22,  # Zn
        1.22,  # Ga
        1.20,  # Ge
        1.19,  # As
        1.20,  # Se
        1.20,  # Br
        1.16,  # Kr
        2.20,  # Rb
        1.95,  # Sr
        1.90,  # Y
        1.75,  # Zr
        1.64,  # Nb
        1.54,  # Mo
        1.47,  # Tc
        1.46,  # Ru
        1.42,  # Rh
        1.39,  # Pd
        1.45,  # Ag
        1.44,  # Cd
        1.42,  # In
        1.39,  # Sn
        1.39,  # Sb
        1.38,  # Te
        1.39,  # I
        1.40,  # Xe
        2.44,  # Cs
        2.15,  # Ba
        2.07,  # La
        2.04,  # Ce
        2.03,  # Pr
        2.01,  # Nd
        1.99,  # Pm
        1.98,  # Sm
        1.98,  # Eu
        1.96,  # Gd
        1.94,  # Tb
        1.92,  # Dy
        1.92,  # Ho
        1.89,  # Er
        1.90,  # Tm
        1.87,  # Yb
        1.87,  # Lu
        1.75,  # Hf
        1.70,  # Ta
        1.62,  # W
        1.51,  # Re
        1.44,  # Os
        1.41,  # Ir
        1.36,  # Pt
        1.36,  # Au
        1.32,  # Hg
        1.45,  # Tl
        1.46,  # Pb
        1.48,  # Bi
        1.40,  # Po
        1.50,  # At
        1.50,  # Rn
        2.60,  # Fr
        2.21,  # Ra
        2.15,  # Ac
        2.06,  # Th
        2.00,  # Pa
        1.96,  # U
        1.90,  # Np
        1.87,  # Pu
        1.80,  # Am
        1.69,  # Cm
        missing,  # Bk
        missing,  # Cf
        missing,  # Es
        missing,  # Fm
        missing,  # Md
        missing,  # No
        missing,  # Lr
        missing,  # Rf
        missing,  # Db
        missing,  # Sg
        missing,  # Bh
        missing,  # Hs
        missing,  # Mt
        missing,  # Ds
        missing,  # Rg
        missing,  # Cn
        missing,  # Nh
        missing,  # Fl
        missing,  # Mc
        missing,  # Lv
        missing,  # Ts
        missing,  # Og
    ]
)


def com(coordinates, atoms):
    masses = np.array([symbol_to_mass[number_to_symbol[a]] for a in atoms], dtype=float)
    com = np.sum(np.multiply(masses[:, np.newaxis], coordinates)) / np.sum(masses)
    return coordinates - com


def calc_pmoi(coordinates, atoms):
    x, y, z = coordinates.T
    masses = np.array([symbol_to_mass[number_to_symbol[a]] for a in atoms], dtype=float)
    Ixx = np.sum(masses * (y**2 + z**2))
    Iyy = np.sum(masses * (x**2 + z**2))
    Izz = np.sum(masses * (x**2 + y**2))
    Ixy = -np.sum(masses * x * y)
    Iyz = -np.sum(masses * y * z)
    Ixz = -np.sum(masses * x * z)
    I0 = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
    I, V = np.linalg.eig(I0)
    idx = np.argsort(-I)
    return V[:, idx].T


def unit_vector(vec):
    return vec / np.linalg.norm(vec)


def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def rotaxis(a, b):
    if np.allclose(a, b):
        return np.array([1, 0, 0])
    c = np.cross(a, b)
    return c / np.linalg.norm(c)


def at_eq(a, b):
    return (a["atomic_number"] == b["atomic_number"]) and (a["degree"] == b["degree"])


def b_eq(a, b):
    return np.isclose(a["distance"], b["distance"], rtol=0.05, atol=0.05)


def center_coordinates(coordinates, atoms):
    coordinates = com(coordinates, atoms)
    for axidx in range(2):
        pmoi = calc_pmoi(coordinates, atoms)
        axis_vector = np.zeros(3)
        axis_vector[axidx] = 1
        theta = angle_between(pmoi[axidx], axis_vector)
        ax = rotaxis(pmoi[axidx], axis_vector)
        rot = scipy.spatial.transform.Rotation.from_rotvec(ax, theta)
        coordinates = rot.apply(coordinates)
    return coordinates


class Molecule:
    def __init__(
        self,
        atoms=None,
        coordinates=None,
        energy=None,
        filename=None,
        name=None,
        lines=None,
        radii=None,
        scale_factor=1.10,
        noh=True,
    ):
        self.scale_factor = scale_factor
        self.radii = radii
        self.name = name
        if filename is not None:
            self.name = splitext(filename)[0]
            self.from_file(filename, noh)
        elif lines is not None:
            self.from_lines(lines, noh)
            self.name = name
        else:
            self.atoms = atoms
            self.coordinates = center_coordinates(coordinates, atoms)
            self.name = name
            if self.radii is None and self.atoms is not None:
                self.set_radii()
            else:
                self.radii = radii
            self.set_am()
            self.set_graph()

    def from_file(self, filename, noh=True):
        f = open(filename, "r")
        V = list()
        atoms = list()
        n_atoms = 0

        # Read the first line to obtain the number of atoms to read
        try:
            n_atoms = int(f.readline().strip())
        except ValueError:
            raise InputError(
                f"Could not obtain the number of atoms in the .xyz file {filename} from first line. Please check format."
            )

        # The title line may contain an energy
        title = f.readline().strip()
        energy = None
        posname = None
        if "scf done:" in title and energy is None:
            try:
                etitle = title.replace(" ", "").split(":")[1].rstrip()
                energy = float(etitle) * ha_to_kcalmol
            except ValueError:
                energy = None
            except AttributeError:
                energy = None
        if "energy:" in title and energy is None:
            try:
                etitle = title.split(":")[1].split(" ")[1].rstrip()
                energy = float(etitle) * ha_to_kcalmol
            except ValueError:
                energy = None
            except AttributeError:
                energy = None
        if "Energy:" in title and energy is None:
            try:
                etitle = title.split(":")[1].split(" ")[1].rstrip()
                energy = float(etitle)  # * ha_to_kcalmol
            except ValueError:
                energy = None
            except AttributeError:
                energy = None
        if "Eopt" in title and energy is None:
            try:
                etitle = title.split("Eopt")[-1].rstrip()
                energy = float(etitle) * ha_to_kcalmol
            except ValueError:
                energy = None
            except AttributeError:
                energy = None
            try:
                posname = title.split(" ")[0].split("/")[-1].split(".")[-1].rstrip()
            except Exception:
                pass
        if energy is None:
            try:
                energy = float(title) * ha_to_kcalmol
            except ValueError:
                energy = None
            except AttributeError:
                energy = None

        # Use the number of atoms to not read beyond the end of a file
        for lines_read, line in enumerate(f):
            if lines_read == n_atoms:
                break

            atom = re.findall(r"[a-zA-Z]+", line)[0]
            atom = symbol_to_number[atom]

            numbers = re.findall(r"[-]?\d+\.\d*(?:[Ee][-\+]\d+)?", line)
            numbers = [float(number) for number in numbers]

            # The numbers are not valid unless we obtain exacly three
            if len(numbers) == 3:
                V.append(np.array(numbers))
                atoms.append(atom)
            else:
                raise InputError(
                    "Reading the .xyz file failed in line {0}. Please check the format.".format(
                        lines_read + 2
                    )
                )

        f.close()
        atoms = np.array(atoms, dtype=int)
        V = center_coordinates(np.array(V), atoms)
        self.title = title
        self.coordinates_with_h = V
        self.atoms_with_h = atoms
        if noh:
            self.atoms = atoms[np.where(atoms > 1)]
            self.coordinates = V[np.where(atoms > 1)]
        else:
            self.atoms = atoms
            self.coordinates = V
        self.energy = energy
        if self.name is None and self.name != posname:
            self.name = posname
        if self.radii is None and self.atoms is not None:
            self.set_radii()
        self.set_am()
        self.set_graph()

    def from_lines(self, lines, noh=True):
        V = list()
        atoms = list()
        n_atoms = 0
        lines_iter = iter(lines)

        # Read the first line to obtain the number of atoms to read
        try:
            n_atoms = int(next(lines_iter).strip())
        except ValueError:
            raise InputError(
                f"Could not obtain the number of atoms in the .xyz file line {lines[0]}. Please check format."
            )

        # The title line may contain an energy
        title = next(lines_iter).strip()
        energy = None
        posname = None
        if "energy:" in title and energy is None:
            try:
                etitle = title.split(":")[1].split(" ")[1].rstrip()
                energy = float(etitle) * ha_to_kcalmol
            except ValueError:
                energy = None
            except AttributeError:
                energy = None
        if "Energy:" in title and energy is None:
            try:
                etitle = title.split(":")[1].split(" ")[1].rstrip()
                energy = float(etitle)  # * ha_to_kcalmol
            except ValueError:
                energy = None
            except AttributeError:
                energy = None
        if "Eopt" in title and energy is None:
            try:
                etitle = title.split("Eopt")[-1].rstrip()
                energy = float(etitle) * ha_to_kcalmol
            except ValueError:
                energy = None
            except AttributeError:
                energy = None
            try:
                posname = title.split(" ")[0].split("/")[-1].split(".")[-1].rstrip()
            except Exception:
                pass
        if energy is None:
            try:
                energy = float(title) * ha_to_kcalmol
            except ValueError:
                energy = None
            except AttributeError:
                energy = None

        # Use the number of atoms to not read beyond the end of a file
        for lines_read, line in enumerate(lines_iter):
            if lines_read == n_atoms:
                break

            atom = re.findall(r"[a-zA-Z]+", line)[0]
            atom = symbol_to_number[atom]

            numbers = re.findall(r"[-]?\d+\.\d*(?:[Ee][-\+]\d+)?", line)
            numbers = [float(number) for number in numbers]

            # The numbers are not valid unless we obtain exacly three
            if len(numbers) == 3:
                V.append(np.array(numbers))
                atoms.append(atom)
            else:
                raise InputError(
                    f"Understanding the chunk of lines failed at line:\n {line}\nPlease check the format."
                )

        atoms = np.array(atoms, dtype=int)
        V = center_coordinates(np.array(V), atoms)
        self.title = title
        self.coordinates_with_h = V
        self.atoms_with_h = atoms
        if noh:
            self.atoms = atoms[np.where(atoms > 1)]
            self.coordinates = V[np.where(atoms > 1)]
        else:
            self.atoms = atoms
            self.coordinates = V
        self.energy = energy
        if self.name is None and self.name != posname:
            self.name = posname
        if self.radii is None and self.atoms is not None:
            self.set_radii()
        self.set_am()
        self.set_graph()

    def set_radii(self):
        radii = np.array([covalent_radii[i] for i in self.atoms], dtype=float)
        self.radii = radii

    def set_am(self):
        n = len(self.atoms)
        am = np.zeros((n, n), dtype=float)
        row, col = np.triu_indices(n, 1)
        dm = scipy.spatial.distance.pdist(self.coordinates)
        rm = self.scale_factor * scipy.spatial.distance.pdist(
            self.radii.reshape(-1, 1), metric=lambda x, y: x + y
        )
        am[row, col] = am[col, row] = dm - rm
        self.am = am < 0

    def set_graph(self):
        G = nx.from_numpy_array(self.am, create_using=nx.Graph)
        an_dict = {i: self.atoms[i] for i in range(len(self.atoms))}
        coord_dict = {i: self.coordinates[i] for i in range(len(self.atoms))}
        nx.set_node_attributes(G, an_dict, "atomic_number")
        nx.set_node_attributes(G, coord_dict, "coordinates")
        ds = np.zeros((len(G.edges())))
        cs = np.zeros_like(ds)
        idxs_1 = []
        idxs_2 = []
        for i, edge in enumerate(G.edges()):
            ds[i] = np.linalg.norm(
                self.coordinates[edge[0]] - self.coordinates[edge[1]]
            )
            cs[i] = self.atoms[edge[0]] * self.atoms[edge[1]] / ds[i] ** 2
            idxs_1.append(edge[0])
            idxs_2.append(edge[0])
        b_dict = nx.edge_betweenness_centrality(G, normalized=False)
        d_dict = {edge: d for edge, d in zip(b_dict.keys(), ds)}
        c_dict = {edge: c for edge, c in zip(b_dict.keys(), cs)}
        i1_dict = {edge: i for edge, i in zip(b_dict.keys(), idxs_1)}
        i2_dict = {edge: i for edge, i in zip(b_dict.keys(), idxs_2)}
        degree_dict = {i: val for (i, val) in G.degree()}
        nx.set_node_attributes(G, degree_dict, "degree")
        nx.set_edge_attributes(G, b_dict, "betweenness")
        nx.set_edge_attributes(G, d_dict, "distance")
        nx.set_edge_attributes(G, c_dict, "coulomb_term")
        nx.set_edge_attributes(G, i1_dict, "idx_1")
        nx.set_edge_attributes(G, i2_dict, "idx_2")
        self.graph = G

    def write(self, rootname="output"):
        filename = f"{rootname}.xyz"
        f = open(filename, "w+")
        with np.printoptions(
            suppress=True, precision=6, formatter={"float_kind": "{:0.6f}".format}
        ):
            print(f"{len(self.atoms_with_h)}", file=f)
            if self.energy is not None:
                printable_energy = np.round(self.energy * kcalmol_to_ha, decimals=6)
                print(f"{printable_energy}", file=f)
            else:
                print(f"{self.title}", file=f)
            for i, atom in enumerate(self.atoms_with_h):
                print(
                    f"{number_to_symbol[atom]:2}    {np.round(self.coordinates_with_h[i][0],decimals=6): }    {np.round(self.coordinates_with_h[i][1],decimals=6): }    {np.round(self.coordinates_with_h[i][2],decimals=6): }",
                    file=f,
                )
        f.close()

    def update(self, atoms, coordinates):
        self.atoms = atoms
        self.coordinates = coordinates
        self.set_radii()
        self.set_am()
        self.set_graph()

    def update_with_h(self, atoms_with_h, coordinates_with_h):
        self.atoms_with_h = atoms_with_h
        self.coordinates_with_h = center_coordinates(coordinates_with_h, atoms_with_h)
        self.update(
            atoms_with_h[np.where(atoms_with_h > 1)],
            coordinates_with_h[np.where(atoms_with_h > 1)],
        )

    def __iter__(self):
        for value in [
            self.energy,
            self.coordinates,
            self.atoms,
            self.radii,
            self.am,
            self.graph,
        ]:
            yield value


def align_mols(mols, mask):
    n = len(mols)
    pos_ibj = []
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            GM = nx.algorithms.isomorphism.GraphMatcher(
                mols[i].graph, mols[j].graph, at_eq, b_eq
            )
            iviews = [
                (
                    np.array(list(im.keys()), dtype=int),
                    np.array(list(im.values()), dtype=int),
                )
                for im in GM.isomorphisms_iter()
            ]
            min_res = np.inf
            for k, (ibi, ibj) in enumerate(iviews):
                print(f"Trying view {k}...")
                ibi, ibj = ibi[ibi.argsort()], ibj[ibi.argsort()]
                print(mols[i].atoms)
                print(mols[j].atoms[ibj])
                res, rotated_coordinates = kabsch_rmsd(
                    mols[i].coordinates, mols[j].coordinates[ibj]
                )
                print(f"Current RMSD is {res}")
                if res < min_res:
                    print(f"Updating best RMSD...")
                    save_ibj = ibj
                    min_res = res
                    save_rotated_coordinates = rotated_coordinates
                if np.isclose(0, min_res):
                    print(f"Current RMSD is {min_res}, we break out!")
                    break
            pos_ibj.append(save_ibj)
            assert all(mols[i].atoms == mols[j].atoms[save_ibj])
            mols[j].update_with_h(mols[j].atoms[save_ibj], save_rotated_coordinates)


def rmsd_matrix(mols, sort=False, truesort=False, normalize=True):
    """
    Compute pairwise RMSD matrix between all molecules in mols.

    Parameters
    ----------
    mols : list of N molecule objects

    Returns
    -------
    M : array
        (N,N) matrix
    maxval : maximum pairwise RMSD
    """
    n = len(mols)
    M = np.zeros((n, n))
    pos_ibj = []
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            if sort:
                GM = nx.algorithms.isomorphism.GraphMatcher(
                    mols[i].graph, mols[j].graph, at_eq, b_eq
                )
                iviews = [
                    (
                        np.array(list(im.keys()), dtype=int),
                        np.array(list(im.values()), dtype=int),
                    )
                    for im in GM.isomorphisms_iter()
                ]
                min_res = np.inf
                for k, (ibi, ibj) in enumerate(iviews):
                    ibi, ibj = ibi[ibi.argsort()], ibj[ibi.argsort()]
                    res, rotated_coordinates = kabsch_rmsd(
                        mols[i].coordinates, mols[j].coordinates[ibj]
                    )
                    if res < min_res:
                        save_ibj = ibj
                        min_res = res
                        save_rotated_coordinates = rotated_coordinates
                    if np.isclose(0, min_res):
                        break
                pos_ibj.append(save_ibj)
                assert all(mols[i].atoms == mols[j].atoms[save_ibj])
                mols[j].update(mols[j].atoms[save_ibj], save_rotated_coordinates)
                M[i, j] = M[j, i] = min_res
            if not sort and len(pos_ibj) > 0:
                for k, ibj in enumerate(pos_ibj):
                    res, rotated_coordinates = kabsch_rmsd(
                        mols[i].coordinates, mols[j].coordinates[ibj]
                    )
                    if res < min_res:
                        save_ibj = ibj
                        min_res = res
                        save_rotated_coordinates = rotated_coordinates
                    if np.isclose(0, min_res):
                        break
                assert all(mols[i].atoms == mols[j].atoms[save_ibj])
                mols[j].update(mols[j].atoms[save_ibj], save_rotated_coordinates)
                M[i, j] = M[j, i] = min_res
            if not sort and len(pos_ibj) == 0:
                res, rotated_coordinates = kabsch_rmsd(
                    mols[i].coordinates, mols[j].coordinates
                )
                assert all(mols[i].atoms == mols[j].atoms)
                mols[j].update(mols[j].atoms, rotated_coordinates)
                M[i, j] = M[j, i] = res
        if not truesort:
            sort = False
    maxval = np.max(M)
    if normalize:
        M = np.abs(M) / maxval
    return M, maxval


def kabsch_rmsd(P, Q):
    """
    Rotate matrix P unto Q using Kabsch algorithm and calculate the RMSD.

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        root-mean squared deviation
    Q : array
        (N,D) matrix, where N is points and D is dimension,
        rotated
    """
    Q = kabsch_rotate(Q, P)
    return rmsd(P, Q), Q


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        Root-mean-square deviation

    """
    return np.sum((V - W) ** 2) / len(V)


def kabsch_rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm.

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    P : array
        (N,D) matrix, where N is points and D is dimension,
        rotated

    """
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    U : matrix
        Rotation matrix (D,D)

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


def quaternion_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD

    based on doi:10.1016/1049-9660(91)90036-O

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    P : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
    """
    rot = quaternion_rotate(P, Q)
    P = np.dot(P, rot)
    return rmsd(P, Q)


def quaternion_transform(r):
    """
    Get optimal rotation
    note: translation will be zero when the centroids of each molecule are the
    same
    """
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot = Wt_r.dot(Q_r)[:3, :3]
    return rot


def makeW(r1, r2, r3, r4=0):
    """
    matrix involved in quaternion rotation
    """
    W = np.asarray(
        [[r4, r3, -r2, r1], [-r3, r4, r1, r2], [r2, -r1, r4, r3], [-r1, -r2, -r3, r4]]
    )
    return W


def makeQ(r1, r2, r3, r4=0):
    """
    matrix involved in quaternion rotation
    """
    Q = np.asarray(
        [[r4, -r3, r2, r1], [r3, r4, -r1, r2], [-r2, r1, r4, r3], [-r1, -r2, -r3, r4]]
    )
    return Q


def quaternion_rotate(X, Y):
    """
    Calculate the rotation

    Parameters
    ----------
    X : array
        (N,D) matrix, where N is points and D is dimension.
    Y: array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rot : matrix
        Rotation matrix (D,D)
    """
    N = X.shape[0]
    W = np.asarray([makeW(*Y[k]) for k in range(N)])
    Q = np.asarray([makeQ(*X[k]) for k in range(N)])
    Qt_dot_W = np.asarray([np.dot(Q[k].T, W[k]) for k in range(N)])
    W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
    A = np.sum(Qt_dot_W, axis=0)
    eigen = np.linalg.eigh(A)
    r = eigen[1][:, eigen[0].argmax()]
    rot = quaternion_transform(r)
    return rot


def centroid(X):
    """
    Calculate the centroid from a vectorset X.

    https://en.wikipedia.org/wiki/Centroid
    Centroid is the mean position of all the points in all of the coordinate
    directions.

    C = sum(X)/len(X)

    Parameters
    ----------
    X : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    C : float
        centeroid

    """
    C = X.mean(axis=0)
    return C


def reorder_hungarian(
    p_atoms,
    q_atoms,
    p_coord,
    q_coord,
):
    """
    Re-orders the input atom list and xyz coordinates using the Hungarian
    method (using optimized column results)
    Parameters
    ----------
    p_atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    q_atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    p_coord : array
        (N,D) matrix, where N is points and D is dimension
    q_coord : array
        (N,D) matrix, where N is points and D is dimension
    Returns
    -------
    view_reorder : array
             (N,1) matrix, reordered indexes of atom alignment based on the
             coordinates of the atoms
    """

    # hview = reorder_hungarian(
    #    mols[i].atoms,
    #    mols[j].atoms,
    #    mols[i].coordinates,
    #    mols[j].coordinates,
    # )
    # mols[j].update(mols[j].atoms[hview], mols[j].coordinates[hview])

    # Find unique atoms
    unique_atoms = np.unique(p_atoms)

    # generate full view from q shape to fill in atom view on the fly
    view_reorder = np.zeros(q_atoms.shape, dtype=int)
    view_reorder -= 1

    for atom in unique_atoms:
        (p_atom_idx,) = np.where(p_atoms == atom)
        (q_atom_idx,) = np.where(q_atoms == atom)

        A_coord = p_coord[p_atom_idx]
        B_coord = q_coord[q_atom_idx]

        view = hungarian(A_coord, B_coord)
        view_reorder[p_atom_idx] = q_atom_idx[view]

    return view_reorder


def hungarian(A, B):
    """
    Hungarian reordering.
    Assume A and B are coordinates for atoms of SAME type only
    """

    distances = cdist(A, B, "euclidean")

    # Perform Hungarian analysis on distance matrix between atoms of 1st
    # structure and trial structure
    indices_a, indices_b = linear_sum_assignment(distances)

    return indices_b


def reorder_distance(
    p_atoms,
    q_atoms,
    p_coord,
    q_coord,
):
    """
    Re-orders the input atom list and xyz coordinates by atom type and then by
    distance of each atom from the centroid.
    Parameters
    ----------
    atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    coord : array
        (N,D) matrix, where N is points and D is dimension
    Returns
    -------
    atoms_reordered : array
        (N,1) matrix, where N is points holding the ordered atoms' names
    coords_reordered : array
        (N,D) matrix, where N is points and D is dimension (rows re-ordered)
    """

    # Find unique atoms
    unique_atoms = np.unique(p_atoms)

    # generate full view from q shape to fill in atom view on the fly
    view_reorder = np.zeros(q_atoms.shape, dtype=int)

    for atom in unique_atoms:
        (p_atom_idx,) = np.where(p_atoms == atom)
        (q_atom_idx,) = np.where(q_atoms == atom)

        A_coord = p_coord[p_atom_idx]
        B_coord = q_coord[q_atom_idx]

        # Calculate distance from each atom to centroid
        A_norms = np.linalg.norm(A_coord, axis=1)
        B_norms = np.linalg.norm(B_coord, axis=1)

        reorder_indices_A = np.argsort(A_norms)
        reorder_indices_B = np.argsort(B_norms)

        # Project the order of P onto Q
        translator = np.argsort(reorder_indices_A)
        view = reorder_indices_B[translator]
        view_reorder[p_atom_idx] = q_atom_idx[view]

    return view_reorder


if __name__ == "__main__":
    # We read files including hydrogens into the coordinates
    A = Molecule(filename=sys.argv[1], noh=False)
    B = Molecule(filename=sys.argv[2], noh=False)

    # Align the molecule objects
    align_mols(mols=[A, B], mask=np.ones((A.coordinates.shape[0])))

    A.write("A")
    B.write("B")
