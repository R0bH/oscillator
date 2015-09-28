"""
A selection of functions designed to operate on lammps files.
"""
import numpy as np
import sys


def get_data_file(input_file_):
    """
    Get the file path to the data file from laamps input.
    Slightly hacked as can't get filter to work...
    """
    with open(input_file_) as f:
        content = f.readlines()
        for i in range(len(content)):
            split = str.split(content[i])
            if len(split) > 0:
                if split[0] == "read_data":
                    data_file = split[1]
        try:
            data_file
        except NameError:
            sys.exit("Error: Can't find data file containing atomic information.")
    return data_file

def get_units(input_file_):
    """
    Get units from lammps input file.
    Currently only "metal" and "Lennard-Jones" units may be used
    Slightly hacked as can't get filter to work...
    """

    with open(input_file_) as input_file:
        read_input = input_file.readlines()
    for i in xrange(len(read_input)):
        split = str.split(read_input[i])
        if len(split) > 0:
            if split[0] == "units":
                units = split[1]
    try:
        units
    except NameError:
        sys.exit("Error: Units not defined in lammps input file.")

    return units

def get_atom_types(input_file_):
    """
    Get atomic species from Lammps input file.
    This function does presume a specific format.
    If atom_types is nonsense compare split index with the input file.
    """
    atom_types = []
    with open(input_file_) as input_file:
        read_input = input_file.readlines()
    for i in xrange(len(read_input)):
        split = str.split(read_input[i])
        if len(split) > 0:
            if split[0] == "pair_coeff":
                for i in xrange(len(split)):
                    if i > 3:
                        atom_types.append(split[i])
    try:
        atom_types
    except NameError:
        sys.exit("Error: Atom types cannot be defined from potential in lammps input file.")
    return atom_types


def get_potential(input_file_):
    """
    Get potential style  and file name from Lammps input file.
    This function does presume a specific format.
    If pair_style is nonsense compare split index with the input file.
    """
    pair_style = []
    with open(input_file_) as input_file:
        read_input = input_file.readlines()
    for i in xrange(len(read_input)):
        split = str.split(read_input[i])
        if len(split) > 0:
            if split[0] == "pair_style":
                pair_style.append(split[1])
            if split[0] == "pair_coeff":
                pair_style.append(split[3])
    try:
        pair_style
    except NameError:
        sys.exit("Error: Pair potential style not defined in lammps input file.")
    return pair_style


def get_radius(lmp_):
    """
    Calculate radius based of box length.
    """
    # Extract volume of simulation cell.
    lmp_.command("variable v equal vol")
    volume = lmp_.extract_variable("v", None, 0)
    # Set move radius to side length.
    # Note: Assumes that the simulation box is cubic.
    return np.power(volume, 1./3.)


def adaptive_radius(accept_, radius_, radius_max_):
    """
    Adaptively adjust the step sized used when generating configurations.
    This will attempt to keep the acceptance ratio at 50%.
    All these numbers are somewhat arbitrary.
    """
    if accept_ < 0.5:
        radius_ *= 0.8
    if accept_ > 0.5:
        radius_ /= 0.8
    if radius_ > radius_max_:
        radius_ = radius_max_
    return radius_


def pos_to_array(lmp_):
    """
    Convert lammps position object into a numpy array.
    """
    num_atoms = lmp_.extract_global("natoms", 0)
    positions = lmp_.extract_atom("x", 3)
    array = np.zeros(3*num_atoms)
    for i in xrange(num_atoms):
        array[3*i] = positions[i][0] 
        array[3*i+1] = positions[i][1]
        array[3*i+2] = positions [i][2]
    return array


def array_to_pos(lmp_,array_):
    """
    Set lammps position object to be the value of a numpy array
    """
    num_atoms = lmp_.extract_global("natoms", 0)
    positions = lmp_.extract_atom("x", 3)
    try:
        for i in xrange(num_atoms):
            positions[i][0] = array_[3*i] 
            positions[i][1] = array_[3*i+1] 
            positions[i][2] = array_[3*i+2]
    except IndexError:
        sys.exit("Error: Dimensions of array and lammps positions do not match.")
