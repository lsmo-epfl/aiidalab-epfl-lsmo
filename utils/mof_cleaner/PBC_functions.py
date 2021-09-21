import copy
import itertools

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.spatial import distance

from .atomic import (
    COVALENT_RADII,
    MASS,
    METALS,
    actinides,
    alkali,
    alkaline_earth,
    lanthanides,
    main_group,
    metalloids,
    metals,
    noble_gases,
    non_metals,
    organic,
    transition_metals,
)

deg2rad = np.pi / 180.0


def readcif(name):
    # FIXME @SMM: improve the function to more general cases of cif files
    with open(name, "r") as fi:
        EIF = fi.readlines()
        cond = True
        cond2 = False
        atom_props_count = 0
        atomlines = []
        counter = 0
        cell_parameter_boundary = [0.0, 0.0]
        for line in EIF:
            line_stripped = line.strip()
            if (not line) or line_stripped.startswith("#"):
                continue
            line_splitted = line.split()
            if line_stripped.startswith("_cell_length_a"):
                cell_a = float(line_splitted[1])
                cell_parameter_boundary[0] = counter + 1
            elif line_stripped.startswith("_cell_length_b"):
                cell_b = float(line_splitted[1])
            elif line_stripped.startswith("_cell_length_c"):
                cell_c = float(line_splitted[1])
            elif line_stripped.startswith("_cell_angle_alpha"):
                cell_alpha = float(line_splitted[1])
            elif line_stripped.startswith("_cell_angle_beta"):
                cell_beta = float(line_splitted[1])
            elif line_stripped.startswith("_cell_angle_gamma"):
                cell_gamma = float(line_splitted[1])
                cell_parameter_boundary[1] = counter + 1
            if cond2 == True and line_stripped.startswith("loop_"):
                break
            else:
                if line_stripped.startswith("_atom"):
                    atom_props_count += 1
                    if line_stripped == "_atom_site_label":
                        type_index = atom_props_count - 1
                    elif line_stripped == "_atom_site_fract_x":
                        fracx_index = atom_props_count - 1
                    elif line_stripped == "_atom_site_fract_y":
                        fracy_index = atom_props_count - 1
                    elif line_stripped == "_atom_site_fract_z":
                        fracz_index = atom_props_count - 1
                    elif "charge" in line_stripped:
                        charge_index = atom_props_count - 1

                    cond2 = True
                elif cond2 == True:
                    if len(line_splitted) == atom_props_count:
                        atomlines.append(line)

            counter += 1

        positions = []
        numbers = []
        atomtypes = []
        atoms = []
        for cn, at in enumerate(atomlines):
            ln = at.strip().split()
            positions.append(
                [float(ln[fracx_index]), float(ln[fracy_index]), float(ln[fracz_index])]
            )
            ln[type_index] = ln[type_index].strip("_")
            at_type = "".join([i for i in ln[type_index] if not i.isdigit()])
            atomtypes.append(at_type)

        cpar = np.array([cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma])
        positions = np.array(positions)
        return cpar, atomtypes, positions


def compute_image_flag(cell, fcoord1, fcoord2):
    invcell = np.linalg.inv(cell)
    supercells = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    fcoords = fcoord2 + supercells
    coords = np.array([np.dot(j, cell) for j in fcoords])
    coord1 = np.dot(fcoord1, cell)
    dists = distance.cdist([coord1], coords)
    dists = dists[0].tolist()
    image = dists.index(min(dists))
    return supercells[image]


def linker_length(adjmat, anchors):
    rows, cols = np.where(adjmat == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    max_length = 0
    min_length = 1000
    for i, j in itertools.combinations(anchors, 2):
        max_length = max(len(nx.shortest_path(gr, i, j)) - 1, max_length)
        min_length = min(len(nx.shortest_path(gr, i, j)) - 1, min_length)
    return (min_length, max_length)


def slice_mat(mat, atoms):
    return np.array(mat[np.ix_(list(atoms), list(atoms))])


def ligand_detect(cell, cart_coords, adj_mat, anchorlist):
    invcell = np.linalg.inv(cell)
    fcoords = np.dot(cart_coords, invcell)
    connected_components = [0]
    checked = []
    periodic_images = []
    if 0 in anchorlist:
        periodic_images.append(np.array([0, 0, 0]))
    counter = 0
    while len(connected_components) < len(cart_coords):
        current_node = connected_components[counter]
        for j, v in enumerate(adj_mat[current_node]):
            if v == 1 and (j not in checked) and (j not in connected_components):
                image_flag = compute_image_flag(cell, fcoords[current_node], fcoords[j])
                fcoords[j] += image_flag
                connected_components.append(j)
                checked.append(j)
                if j in anchorlist:
                    periodic_images.append(image_flag)
        counter += 1
    return np.array(periodic_images)


def XYZ_connected(cell, cart_coords, adj_mat):
    invcell = np.linalg.inv(cell)
    fcoords = np.dot(cart_coords, invcell)
    connected_components = [0]
    checked = []
    counter = 0
    while len(connected_components) < len(cart_coords):
        current_node = connected_components[counter]
        for j, v in enumerate(adj_mat[current_node]):
            if v == 1 and j not in checked and j not in connected_components:
                fcoords[j] += compute_image_flag(
                    cell, fcoords[current_node], fcoords[j]
                )
                connected_components.append(j)
                checked.append(j)
        counter += 1
    return fcoords


def writeCIF(filename, cell, atoms_types, fcoords, indices2write=None):
    if not indices2write:
        indices2write = np.arange(len(fcoords))
    [cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma] = make_cellparams(cell)
    with open(filename, "w") as f_cif:
        f_cif.write("# This file is generated after removing free solvents.\n")
        f_cif.write("data_I\n")
        f_cif.write("_cell_length_a %8.04f\n" % (cell_a))
        f_cif.write("_cell_length_b %8.04f\n" % (cell_b))
        f_cif.write("_cell_length_c %8.04f\n" % (cell_c))
        f_cif.write("_cell_angle_alpha %4.02f\n" % (cell_alpha))
        f_cif.write("_cell_angle_beta  %4.02f\n" % (cell_beta))
        f_cif.write("_cell_angle_gamma %4.02f\n" % (cell_gamma))
        f_cif.write("_space_group_name_H-M_alt      'P 1'\n\n\n")
        f_cif.write("loop_\n_space_group_symop_operation_xyz\n  'x, y, z' \n\n")

        f_cif.write("loop_\n")
        f_cif.write("_atom_site_label\n")
        f_cif.write("_atom_site_fract_x\n")
        f_cif.write("_atom_site_fract_y\n")
        f_cif.write("_atom_site_fract_z\n")
        f_cif.write("_atom_site_type_symbol\n")
        for i in indices2write:
            f_cif.write(
                "%-5s %10.4f %10.4f %10.4f %2s\n"
                % (
                    atoms_types[i],
                    fcoords[i][0],
                    fcoords[i][1],
                    fcoords[i][2],
                    atoms_types[i],
                )
            )


def writeXYZfcoords(filename, atoms, cell, fcoords):
    with open(filename, "w") as fo:
        fo.write("%i\n\n" % len(atoms))
        for i, fcoord in enumerate(fcoords):
            cart_coord = np.dot(fcoord, cell)
            s = "%10.2f %10.2f %10.2f" % (cart_coord[0], cart_coord[1], cart_coord[2])
            fo.write("%s %s\n" % (atoms[i], s))


def writeXYZandGraph(filename, atoms, cell, fcoords, molgraph):
    with open(filename, "w") as fo:
        fo.write("%i\n\n" % len(atoms))
        for i, fcoord in enumerate(fcoords):
            cart_coord = np.dot(fcoord, cell)
            s = "%10.2f %10.2f %10.2f" % (cart_coord[0], cart_coord[1], cart_coord[2])
            fo.write("%s %s\n" % (atoms[i], s))
    tmpstr = ",".join([at for at in atoms])
    np.savetxt(filename[:-4] + ".net", molgraph, fmt="%i", delimiter=",", header=tmpstr)


def writeXYZcoords(filename, atoms, coords):
    with open(filename, "w") as fo:
        fo.write("%i\n\n" % len(atoms))
        for i, cart_coord in enumerate(coords):
            s = "%10.2f %10.2f %10.2f" % (cart_coord[0], cart_coord[1], cart_coord[2])
            fo.write("%s %s\n" % (atoms[i], s))
    return


def writeXYZcoords_withcomment(filename, atoms, coords, comment):
    with open(filename, "w") as fo:
        fo.write("%i\n" % len(atoms))
        fo.write("%s\n" % comment)
        for i, cart_coord in enumerate(coords):
            s = "%10.2f %10.2f %10.2f" % (cart_coord[0], cart_coord[1], cart_coord[2])
            fo.write("%s %s\n" % (atoms[i], s))
    return


def write2file(pt, fn, st):
    with open(pt + fn, "a") as fo:
        fo.write(st)


def min_img_distance(coords1, coords2, cell):
    invcell = np.linalg.inv(cell)
    one = np.dot(coords1, invcell) % 1
    two = np.dot(coords2, invcell) % 1
    three = np.around(one - two)
    four = np.dot(one - two - three, cell)
    return np.linalg.norm(four)


def fractional2cart(fcoord, cell):
    return np.dot(fcoord, cell)


def min_img_distance2(coord1, coord2, cell):
    invcell = np.linalg.inv(cell)
    supercells = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    fcoords = np.dot(coord2, invcell) + supercells
    coords = np.array([np.dot(j, cell) for j in fcoords])
    dists = distance.cdist([coord1], coords)
    return np.amin(dists)


def frac_coord(coord, cell):
    invcell = np.linalg.inv(cell)
    return np.dot(coord, invcell)


def compute_distance_matrix(cell, cart_coords):
    distance_matrix = np.zeros([len(cart_coords), len(cart_coords)])
    for i in range(len(cart_coords)):
        for j in range(i + 1, len(cart_coords)):
            d = min_img_distance(cart_coords[i], cart_coords[j], cell)
            distance_matrix[i, j] = d
            distance_matrix[j, i] = d

    return distance_matrix


def compute_distance_matrix2(cell, cart_coords):
    distance_matrix = np.zeros([len(cart_coords), len(cart_coords)])
    for i in range(len(cart_coords)):
        for j in range(i + 1, len(cart_coords)):
            d = min_img_distance2(cart_coords[i], cart_coords[j], cell)
            distance_matrix[i, j] = d
            distance_matrix[j, i] = d

    return distance_matrix


def make_graph_from_nodes_edges(nodes, edges, attribs):
    gr = nx.Graph()
    [gr.add_node(n, atomicNum=at) for n, at in zip(nodes, attribs)]
    # gr.add_nodes_from(nodes)
    gr.add_edges_from(edges)
    return gr


def make_cellparams(cellv):
    A = cellv[0]
    B = cellv[1]
    C = cellv[2]
    cell_a = round(np.sqrt(A[0] ** 2 + A[1] ** 2 + A[2] ** 2), 3)
    cell_b = round(np.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2), 3)
    cell_c = round(np.sqrt(C[0] ** 2 + C[1] ** 2 + C[2] ** 2), 3)
    cell_gamma = round(np.arccos(B[0] / cell_b) / deg2rad, 3)
    cell_beta = round(np.arccos(C[0] / cell_c) / deg2rad, 3)
    cell_alpha = round(
        np.arccos(
            C[1] / cell_c * np.sin(cell_gamma * deg2rad)
            + np.cos(cell_beta * deg2rad) * np.cos(cell_gamma * deg2rad)
        )
        / deg2rad,
        3,
    )
    return [cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma]


def mkcell(params):
    """Update the cell representation to match the parameters."""
    a_mag, b_mag, c_mag = params[:3]
    alpha, beta, gamma = [x * deg2rad for x in params[3:]]
    a_vec = np.array([a_mag, 0.0, 0.0])
    b_vec = np.array([b_mag * np.cos(gamma), b_mag * np.sin(gamma), 0.0])
    c_x = c_mag * np.cos(beta)
    c_y = c_mag * (np.cos(alpha) - np.cos(gamma) * np.cos(beta)) / np.sin(gamma)
    c_vec = np.array([c_x, c_y, (c_mag ** 2 - c_x ** 2 - c_y ** 2) ** 0.5])
    return np.array([a_vec, b_vec, c_vec])


def make_supercell(cell, atoms, fcoords, exp_coeff):
    supercell = np.multiply(cell.T, exp_coeff).T
    superatoms = []
    superfcoords = []
    for i in range(exp_coeff[0]):
        for j in range(exp_coeff[1]):
            for k in range(exp_coeff[2]):
                for na, atom in enumerate(atoms):
                    cpatom = {}
                    fc = fcoords[na]
                    fx = fc[0] / exp_coeff[0] + float(i) / exp_coeff[0]
                    fy = fc[1] / exp_coeff[1] + float(j) / exp_coeff[1]
                    fz = fc[2] / exp_coeff[2] + float(k) / exp_coeff[2]
                    superfcoords.append([fx, fy, fz])
                    superatoms.append(atom)
    superfcoords = np.array(superfcoords)
    return supercell, superatoms, superfcoords


def compute_overlap_matrix(distance_mat, allatomtypes, tolerance=1.0):
    """
    Find atomic overlap based on pairwise distance and Covalent radii.
    Current criterion: if dist < min (CovR_1,CovR_2) -> overlap
    """
    overlap_matrix = np.zeros(distance_mat.shape)
    for i, e1 in enumerate(allatomtypes[:-1]):
        for j, e2 in enumerate(allatomtypes[i + 1 :]):
            dist = distance_mat[i, i + j + 1]
            # check for atomic overlap:
            if dist < tolerance * min(COVALENT_RADII[e1], COVALENT_RADII[e2]):
                overlap_matrix[i, i + j + 1] = 1
                overlap_matrix[i + j + 1, i] = 1
    return sparse.csr_matrix(overlap_matrix)


def compute_adj_matrix(distance_mat, allatomtypes):
    adj_matrix = np.zeros(distance_mat.shape)
    for i, e1 in enumerate(allatomtypes[:-1]):
        for j, e2 in enumerate(allatomtypes[i + 1 :]):
            elements = set([e1, e2])
            if elements < metals:  # FIXME no metal-metal bond allowed
                continue
            rad = COVALENT_RADII[e1] + COVALENT_RADII[e2]
            dist = distance_mat[i, i + j + 1]
            # check for atomic overlap:
            # if dist < min(COVALENT_RADII[e1] , COVALENT_RADII[e2]):
            #     print("atomic overlap!")
            #     raise NotImplementedError
            tempsf = 0.9
            # probably a better way to fix these kinds of issues..
            if (set("F") < elements) and (elements & metals):
                tempsf = 0.8
            if (set("C") < elements) and (elements & metals):
                tempsf = 0.95
            if (
                (set("H") < elements)
                and (elements & metals)
                and (not elements & alkali)
            ):
                tempsf = 0.75

            if (set("O") < elements) and (elements & metals):
                tempsf = 0.85
            if (set("N") < elements) and (elements & metals):
                tempsf = 0.82
            # fix for water particle recognition.
            if set(["O", "H"]) <= elements:
                tempsf = 0.8
            # very specific fix for Michelle's amine appended MOF
            if set(["N", "H"]) <= elements:
                tempsf = 0.67
            if set(["Mg", "N"]) <= elements:
                tempsf = 0.80
            if set(["C", "H"]) <= elements:
                tempsf = 0.80
            if set(["K"]) <= elements:
                tempsf = 0.95
            if lanthanides & elements:
                tempsf = 0.95
            if elements == set(["C"]):
                tempsf = 0.85
            if dist * tempsf < rad:  # and not (alkali & elements):
                adj_matrix[i, i + j + 1] = 1
                adj_matrix[i + j + 1, i] = 1
    return sparse.csr_matrix(adj_matrix)


def get_closed_subgraph(linkers, SBUlist, adj_matrix):
    ###############################################################################
    # This part separates the linkers into their respective subgraphs             #
    # First element is the things you want to find subgraphs of.                  #
    # If this is the linkers, you input that as the first.                        #
    # If you input the SBU as the first, then you get the subgraphs of the SBU.   #
    # The second element tells you what part of the matrix is NOT what you want.  #
    # If we want subgraphs of linkers, we want to exclude the SBU.                #
    ###############################################################################

    linkers_sub = linkers.copy()
    linker_list = []
    linker_subgraphlist = []
    counter = 0
    while len(linkers_sub) > 0:
        counter += 1
        if counter > 5000:
            break
        start_idx = list(linkers_sub)[0]
        current_linker_list = set([start_idx])
        checked_list = set()
        while len(checked_list) <= len(current_linker_list):
            loop_over = np.nonzero(adj_matrix[start_idx])[1]
            current_linker_list.update(np.nonzero(adj_matrix[start_idx])[1])
            current_linker_list = current_linker_list - SBUlist
            checked_list.add(start_idx)
            for val in loop_over:
                if not val in SBUlist:
                    current_linker_list.update(np.nonzero(adj_matrix[val])[1])
            left_to_check = current_linker_list - checked_list - SBUlist
            if len(left_to_check) == 0:
                break
            else:
                start_idx = list(left_to_check)[0]
        current_linker_list = current_linker_list - SBUlist
        linkers_sub = linkers_sub - current_linker_list
        ####### We want to return both the linker itself as well as the subgraph corresponding to it.
        linker_list.append(list(current_linker_list))
        linker_subgraphlist.append(
            adj_matrix[np.ix_(list(current_linker_list), list(current_linker_list))]
        )
    return linker_list, linker_subgraphlist


def include_extra_shells(SBUlists, subgraphlists, molcif, adjmat):
    SBUs = []
    subgraphs = []
    for SBU in SBUlists:
        for zero_first_shell in copy.deepcopy(SBU):
            for val in molcif.getBondedAtomsSmart(zero_first_shell):
                SBU.append(val)
        SBUset = set(SBU)
        SBUs.append(list(SBUset))
        subgraphs.append(adjmat[np.ix_(list(SBUset), list(SBUset))])

    return SBUs, subgraphs


def find_metal(atom_types):
    metallist = []
    for i, at in enumerate(atom_types):
        if at in metals:
            metallist.append(i)

    return metallist
