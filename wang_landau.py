#!/usr/bin/python
"""
Performs Wang-Landau Monte Carlo using Lammps as an energy solver.
"""

import os
import subprocess
import numpy as np
from lammps_commands import get_radius
from atom_commands import atoms
                          


def check_lnf(histogram_, update_lnf_):
    """
    Checks if each histogram region has been visited 1/sqrt(lnf) times.
    Used as a convergence criterion for WL.
    Note: This is the default convergence function.
    """
    min_val = np.amin(histogram_[histogram_ > 0])
    if min_val > 1./np.sqrt(update_lnf_):
        return True
    else:
        return False


def check_flatness(histogram_):
    """
    Checks if a histogram is flat to a given tolerance.
    Used as a convergence criterion for WL.
    This is the original Wang-Landau criterion.
    Note: This is not default must be uncommented to use.
    """
    mean = np.average(histogram_[histogram_ > 0])
    min_val = np.amin(histogram_[histogram_ > 0])
    if min_val > 0.8*mean:
        return True
    else:
        return False


def wang_landau(atom_system_, lmp_, update_lnf_, density_of_states_=None,
                e_lower_=-65., e_higher_=22.):
    """
    Performs a Wang - Landau simulation on a lammps instance (lmp).
    The WL simulation is carried out with update value update_lnf_.
    The simulation can be started from a previous guess of the DOS.
    """
    # Define commands for optimisation purposes.
    gen_rand = np.random.rand
    pick_atom = atom_system_.get_random_atom
    move = atom_system_.move_atom
    check_converge = check_lnf

    # Open path for Wang-Landau output.
    root_folder = "wl_output"
    if not os.path.exists(root_folder):
        os.mkdir(root_folder)
    path = root_folder+"/wl_output.txt"
    wl_output = open(path, 'a')

    # Define upper and lower bound for the energies of the system.
    # A range is then created within these bounds.
    # Note: Only do this if the DOS has not already been defined.
    num_bins = 500
    spacing = (e_higher_ - e_lower_)/num_bins
    e_range = np.arange(e_lower_, e_higher_, spacing)
    # Define histogram based on the above defined range.
    histogram = np.zeros(len(e_range))
    # Define minimum number of iterations.
    # Calculated from the number of bins multiplied by min visits to each bin.
    min_its = num_bins/(np.sqrt(update_lnf_))
    if density_of_states_ is None:
        print "Note: Initialising DOS for lnf = "+str(update_lnf_)
        wl_output.write("Note: Initialising DOS for lnf = " \
                        +str(update_lnf_)+"\n")
        density_of_states_ = np.zeros(len(e_range))
    else:
        print "lnf = "+str(update_lnf_)
        wl_output.write("lnf = "+str(update_lnf_)+"\n")

    # Define total number of WL iterations.
    wl_its = 500000000

    # Run lammps to calculate inital properties of the system.
    lmp_.command("run 0 post no")
    lmp_.command("variable e equal pe")

    # Extract the atomic coordinates.
    # This is an array of pointers.
    # Note: Changing these values will alter the values in lammps.
    atom_system_.positions = lmp_.extract_atom("x", 3)

    # Set move radius to side length.
    # Note: Assumes that the simulation box is cubic.
    radius = get_radius(lmp_)/20.

    atom_system_.generate_configuration(lmp_, radius, lower_bound_=e_lower_,
                           upper_bound_=e_higher_)

    # Calculate and store the total energy of the initial configuration.
    lmp_.command("variable elast equal $e")
    elast = lmp_.extract_variable("elast", None, 0)
    print "Energy of initial configuration = " + str(elast)
    wl_output.write("Energy of initial configuration = " + str(elast)+"\n")
    # Coordinate store for MC.
    positions_store = np.zeros(3)

    # acceptance counter (eventually implement an adaptive step size)
    accept = 0

    lmp_.command("dump id_wl_initial all xyz 50 wl_initial_state.xyz")
    lmp_.command("run 0 post no")
    lmp_.command("undump id_wl_initial")

    # Initialise level
    # This needs reworking as it breaks if initial energy lies outside range.
    old_level = int((elast - e_lower_)/spacing)

    # Wang-Landau Monte Carlo Algorithm.
    for i in xrange(1, wl_its):

        # Select an atom at random.
        atom = pick_atom()
        # Store previous positions.
        positions_store = [atom_system_.positions[atom][0],
                           atom_system_.positions[atom][1],
                           atom_system_.positions[atom][2]]
        # Move selected atom to generate a new configuration.
        move(atom, radius)
        # Run lammps for the new configuration.
        lmp_.command("run 0 post no")
        # Refresh positions from lammps.
        # Required due to book-keeping (neighbour lists etc).
        atom_system_.positions = lmp_.extract_atom("x", 3)
        # Extract the total energy of the new configuration.
        lmp_.command("variable enew equal v_e")
        enew = lmp_.extract_variable("enew", None, 0)
        # Apply WL criterion to see if the move is accepted.
        level = int((elast - e_lower_)/spacing)
        if enew < e_higher_ and enew > e_lower_:
            if gen_rand() < \
               np.exp(density_of_states_[level]-density_of_states_[old_level]):
                # If the move is accepted update stored energy.
                elast = enew
                lmp_.command("variable elast equal v_e")
                old_level = level
                # update acceptance counter
                accept += 1
            else:
                # If the move is rejected revert to original position of atom.
                atom_system_.positions[atom][0] = positions_store[0]
                atom_system_.positions[atom][1] = positions_store[1]
                atom_system_.positions[atom][2] = positions_store[2]
            # Update histogram in energy space recording final configuration.
            histogram[level] += 1
            density_of_states_[level] += update_lnf_
        else:
            # If the move is rejected revert to original position of atom.
            atom_system_.positions[atom][0] = positions_store[0]
            atom_system_.positions[atom][1] = positions_store[1]
            atom_system_.positions[atom][2] = positions_store[2]
        # Check convergence
        if (i > min_its) and (i % 10000 == 0):
            if check_converge(histogram, update_lnf_):
                print "Converged after " + str(i+1)+" iterations."
                wl_output.write("Converged after " + str(i+1)+" iterations.\n")
                print "Acceptance = " + str(float(accept)/float(i))+"\n"
                wl_output.write("Acceptance = " \
                                + str(float(accept)/float(i))+"\n")
                path = root_folder+'/wl_dos_'+str(update_lnf_)+'.txt'
                write_dos = open(path, 'w')
                for i in range(len(histogram)):
                    write_dos.write(str(e_range[i]) +
                                    '    ' + str(density_of_states_[i]) +
                                    '    ' + str(histogram[i]) +
                                    '\n')
                return density_of_states_, True

    return density_of_states_, False
