import re
import sys
from os.path import dirname, splitext
import networkx as nx
import numpy as np
import scipy.spatial
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import sys

ha_to_kcalmol = 627.509474
kcalmol_to_ha = 0.00159360
symbol_to_number = {
    "Em": 0,
    "Vc": 0,
    "Va": 0,
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
        missing,
        0.31,
        0.28,
        1.28,
        0.96,
        0.84,
        0.76,
        0.71,
        0.66,
        0.57,
        0.58,
        1.66,
        1.41,
        1.21,
        1.11,
        1.07,
        1.05,
        1.02,
        1.06,
        2.03,
        1.76,
        1.70,
        1.60,
        1.53,
        1.39,
        1.39,
        1.32,
        1.26,
        1.24,
        1.32,
        1.22,
        1.22,
        1.20,
        1.19,
        1.20,
        1.20,
        1.16,
        2.20,
        1.95,
        1.90,
        1.75,
        1.64,
        1.54,
        1.47,
        1.46,
        1.42,
        1.39,
        1.45,
        1.44,
        1.42,
        1.39,
        1.39,
        1.38,
        1.39,
        1.40,
        2.44,
        2.15,
        2.07,
        2.04,
        2.03,
        2.01,
        1.99,
        1.98,
        1.98,
        1.96,
        1.94,
        1.92,
        1.92,
        1.89,
        1.90,
        1.87,
        1.87,
        1.75,
        1.70,
        1.62,
        1.51,
        1.44,
        1.41,
        1.36,
        1.36,
        1.32,
        1.45,
        1.46,
        1.48,
        1.40,
        1.50,
        1.50,
        2.60,
        2.21,
        2.15,
        2.06,
        2.00,
        1.96,
        1.90,
        1.87,
        1.80,
        1.69,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
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
        try:
            n_atoms = int(f.readline().strip())
        except ValueError:
            raise InputError(
                f"Could not obtain the number of atoms in the .xyz file {filename} from first line. Please check format."
            )
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
                energy = float(etitle)
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
        for lines_read, line in enumerate(f):
            if lines_read == n_atoms:
                break
            atom = re.findall(r"[a-zA-Z]+", line)[0]
            atom = symbol_to_number[atom]
            numbers = re.findall(r"[-]?\d+\.\d*(?:[Ee][-\+]\d+)?", line)
            numbers = [float(number) for number in numbers]
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
        try:
            n_atoms = int(next(lines_iter).strip())
        except ValueError:
            raise InputError(
                f"Could not obtain the number of atoms in the .xyz file line {lines[0]}. Please check format."
            )
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
                energy = float(etitle)
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
        for lines_read, line in enumerate(lines_iter):
            if lines_read == n_atoms:
                break
            atom = re.findall(r"[a-zA-Z]+", line)[0]
            atom = symbol_to_number[atom]
            numbers = re.findall(r"[-]?\d+\.\d*(?:[Ee][-\+]\d+)?", line)
            numbers = [float(number) for number in numbers]
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
    print(f"Aligning {n} molecules...")
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


def kabsch_weighted_fit(
    P,
    Q,
    W=None,
    return_rmsd=False,
):
    rmsd_: float
    R, T, rmsd_ = kabsch_weighted(Q, P, W)
    PNEW = np.dot(P, R.T) + T
    if return_rmsd:
        return (PNEW, rmsd_)
    else:
        return (PNEW, None)


def kabsch_weighted_rmsd(P, Q, W=None):
    U, V, w_rmsd = kabsch_weighted(P, Q, W)
    return w_rmsd, U, V


def rmsd_matrix(mols, sort=False, truesort=False, normalize=True):
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


def kabsch_rmsd(P, Q, W=None):
    if W is not None:
        return kabsch_weighted_rmsd(P, Q, W)
    U = kabsch_rotate(Q, P)
    return rmsd(P, Q), U


def rmsd(V, W):
    diff = V - W
    return np.sqrt((diff * diff).sum() / V.shape[0])


def kabsch_rotate(P, Q):
    U = kabsch(P, Q)
    return U


def kabsch(P, Q):
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    return U


def kabsch_weighted(P, Q, W=None):
    CMP = np.zeros(3)
    CMQ = np.zeros(3)
    C = np.zeros((3, 3))
    if W is None:
        W = np.ones(len(P)) / len(P)
    W = np.array([W, W, W]).T
    iw = 3.0 / W.sum()
    n = len(P)
    for i in range(3):
        for j in range(n):
            for k in range(3):
                C[i, k] += P[j, i] * Q[j, k] * W[j, i]
    CMP = (P * W).sum(axis=0)
    CMQ = (Q * W).sum(axis=0)
    PSQ = (P * P * W).sum() - (CMP * CMP).sum() * iw
    QSQ = (Q * Q * W).sum() - (CMQ * CMQ).sum() * iw
    C = (C - np.outer(CMP, CMQ) * iw) * iw
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    msd = (PSQ + QSQ) * iw - 2.0 * S.sum()
    if msd < 0.0:
        msd = 0.0
    rmsd_ = np.sqrt(msd)
    V = np.zeros(3)
    for i in range(3):
        t = (U[i, :] * CMQ).sum()
        V[i] = CMP[i] - t
    V = V * iw
    return U, V, rmsd_


def quaternion_rmsd(P, Q):
    rot = quaternion_rotate(P, Q)
    P = np.dot(P, rot)
    return rmsd(P, Q)


def quaternion_transform(r):
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot = Wt_r.dot(Q_r)[:3, :3]
    return rot


def makeW(r1, r2, r3, r4=0):
    W = np.asarray(
        [[r4, r3, -r2, r1], [-r3, r4, r1, r2], [r2, -r1, r4, r3], [-r1, -r2, -r3, r4]]
    )
    return W


def makeQ(r1, r2, r3, r4=0):
    Q = np.asarray(
        [[r4, -r3, r2, r1], [r3, r4, -r1, r2], [-r2, r1, r4, r3], [-r1, -r2, -r3, r4]]
    )
    return Q


def quaternion_rotate(X, Y):
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
    C = X.mean(axis=0)
    return C


def reorder_hungarian(
    p_atoms,
    q_atoms,
    p_coord,
    q_coord,
):
    unique_atoms = np.unique(p_atoms)
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
    distances = cdist(A, B, "euclidean")
    indices_a, indices_b = linear_sum_assignment(distances)
    return indices_b


def reorder_distance(
    p_atoms,
    q_atoms,
    p_coord,
    q_coord,
):
    unique_atoms = np.unique(p_atoms)
    view_reorder = np.zeros(q_atoms.shape, dtype=int)
    for atom in unique_atoms:
        (p_atom_idx,) = np.where(p_atoms == atom)
        (q_atom_idx,) = np.where(q_atoms == atom)
        A_coord = p_coord[p_atom_idx]
        B_coord = q_coord[q_atom_idx]
        A_norms = np.linalg.norm(A_coord, axis=1)
        B_norms = np.linalg.norm(B_coord, axis=1)
        reorder_indices_A = np.argsort(A_norms)
        reorder_indices_B = np.argsort(B_norms)
        translator = np.argsort(reorder_indices_A)
        view = reorder_indices_B[translator]
        view_reorder[p_atom_idx] = q_atom_idx[view]
    return view_reorder


# Ignore everything except this bit maybe


def rmsd_masked(P, Q, indices_fixed, indices_sub, translate=True):

    if translate:
        Q = Q - centroid(Q)
        P = P - centroid(P)

    # We set up indices to only consider substrate in first rotation
    W_fixed = np.zeros_like((P[:, 0]))
    indices_fixed = np.array(indices_fixed, dtype=int)
    W_fixed[indices_fixed] = 1

    # We compute RMSD and rotation bits for fixed part (substrate)
    res_fixed, U, V = kabsch_rmsd(P, Q, W_fixed)
    #print(res_fixed)
    # The RMSD between the two fixed parts should be small since the substrate doesnt change much. This check can be relaxed.
    #assert np.isclose(res_fixed, 0, atol=1.1)#0.5)#1e-02)#4)

    # We apply the rotation
    P = np.dot(P-V, U)
    #print(P[indices_fixed])
    #print(Q[indices_fixed])
    # Once rotated, the coordinates of the fixed part should also be almost identical. This check can be relaxed.
    #assert np.allclose(P[indices_fixed], Q[indices_fixed], atol=1.1)#0.5)#1e-02)

    # We can now compute RMSD for the substituent part
    return rmsd(Q[indices_sub], P[indices_sub])


if __name__ == "__main__":
    from AaronTools.geometry import Geometry
    int1 = sys.argv[1]
    int2 = sys.argv[2]

    A = Geometry(int1)
    B = Geometry(int2)
    #print(A)
    P = A.coordinates()
    Q = B.coordinates()

    # We read files including hydrogens into the coordinates
    # Quick tests
    mask_scaffold = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    mask_substituent = [20]#[29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
    #mask_1 = []
    rmsd_1 = rmsd_masked(P, Q, mask_scaffold, mask_substituent)
    print(rmsd_1)
