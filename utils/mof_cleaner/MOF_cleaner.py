import itertools
import sys
from collections import Counter

import networkx as nx
import numpy as np
from scipy import sparse

from .PBC_functions import (
    compute_adj_matrix,
    compute_overlap_matrix,
    find_metal,
    fractional2cart,
    mkcell,
    writeCIF,
)


def clean_MOF_free_solvents(
    cif,
    moleculesize_ratio=0.5,
    framework_min_size=20,
    solvent_max_size=100,
    check_metal=True,
):
    """Find the free solvent molecules in the pores. Takes a ASE structure and return solvent indices."""

    output_str = ""
    if hasattr(cif, "periodic_distance_matrix"):
        distance_mat = cif.periodic_distance_matrix
    else:
        distance_mat = cif.get_all_distances(mic=True)
        cif.periodic_distance_matrix = distance_mat
    allatomtypes = cif.get_chemical_symbols()
    adj_matrix = compute_adj_matrix(distance_mat, allatomtypes)

    # Check number of connected components.
    # if more than 1: it checks if the structure is interpenetrated. Fails if no metal in one of the connected components (identified by the graph).
    # This includes floating solvent molecules.

    n_components, labels_components = sparse.csgraph.connected_components(
        csgraph=adj_matrix, directed=False, return_labels=True
    )
    # Two scenarios for finding solvents:
    #    1. based on the size of component
    #    2. based on metalic components
    # first we find the largest connected components, if it has metal, we assume the structure is a MOF

    components_sizes = dict(Counter(labels_components))
    largest_component = max(components_sizes, key=components_sizes.get)
    largest_component_size = components_sizes[largest_component]

    # Check if the components based on metal.
    if check_metal:
        metal_list = set(find_metal(allatomtypes))
        cif_has_metal = False
        inds_in_comp = [
            i
            for i in range(len(labels_components))
            if labels_components[i] == largest_component
        ]
        if set(inds_in_comp) & metal_list:
            cif_has_metal = True

    main_components = []
    main_components_atoms = []
    solvent_components = []
    solvent_components_atoms = []
    solvent_id = 0
    for comp in range(n_components):
        inds_in_comp = [
            i for i in range(len(labels_components)) if labels_components[i] == comp
        ]
        if (
            components_sizes[comp] <= moleculesize_ratio * largest_component_size
            or components_sizes[comp] < framework_min_size
        ) and (components_sizes[comp] < solvent_max_size):
            sol_atoms = []
            for i in inds_in_comp:
                solvent_components.append(i)
                solvent_components_atoms.append(allatomtypes[i])
                sol_atoms.append(allatomtypes[i])
            solvent_id += 1
        elif check_metal and cif_has_metal:
            if not set(inds_in_comp) & metal_list:
                sol_atoms = []
                for i in inds_in_comp:
                    solvent_components.append(i)
                    solvent_components_atoms.append(allatomtypes[i])
                    sol_atoms.append(allatomtypes[i])
                solvent_id += 1
        else:
            main_components.append(comp)
            main_components_atoms = main_components_atoms + inds_in_comp

    if solvent_id > 0:
        output_str += "Found {} solvents in the pore\n".format(solvent_id)

    return solvent_components, solvent_components_atoms, output_str


def find_atomic_overlap(cif, tolerance=1.0):
    """Find atomic overlap based on covalent radii.

    # TODO: we can return pairs of atoms with overlap such that the user can go through them
    """
    output_str = ""
    if hasattr(cif, "periodic_distance_matrix"):
        distance_mat = cif.periodic_distance_matrix
    else:
        distance_mat = cif.get_all_distances(mic=True)
        cif.periodic_distance_matrix = distance_mat
    allatomtypes = cif.get_chemical_symbols()
    overlap_matrix = compute_overlap_matrix(distance_mat, allatomtypes, tolerance)
    overlap_atoms = []
    for at in list(set(sparse.find(overlap_matrix)[0])):
        overlap_atoms.append(at.item())
    if not overlap_atoms:
        output_str += "no atomic overlap was found\n"
    else:
        output_str += "atomic overlap was found for {} atoms".format(len(overlap_atoms))
    overlap_atoms_types = [allatomtypes[i] for i in overlap_atoms]
    return overlap_atoms, overlap_atoms_types, output_str
