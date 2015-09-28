"""
Commands used to select and move atoms.
"""

import numpy as np
import sys
from lammps_commands import array_to_pos, pos_to_array, get_data_file, get_atom_types


class atoms:
    """ The atomic information """
    def __init__(self, lmp_, input_file_ = None):
        self.num_atoms = lmp_.extract_global("natoms", 0)
        self.num_atom_types = []
        self.atom_types = []
        self.initial_positions = pos_to_array(lmp_)
        self.positions_store = pos_to_array(lmp_)
        self.positions = pos_to_array(lmp_)
        if input_file_ == None:
            print "No input file defined in atoms system"
            exit()
        else:    
            """
            Assign atom types.
            """
            data_file = get_data_file(input_file_)
            print "Obtaining atom types from " + data_file
            with open(data_file) as f:
                content = f.readlines()
            types_array=get_atom_types(input_file_)
            self.num_atom_types = int(str.split(content[3])[0])
            print str(self.num_atom_types) + " found of type:"
            print types_array
            for i in range(self.num_atoms):
                atom_type = int(str.split(content[14+(self.num_atom_types)+i])[1])
                self.atom_types.append(types_array[atom_type-1])
            print "Atom types are as follows:"
            print self.atom_types
            print "\n"

    def get_num_atoms(self):
        return self.num_atoms

    def set_num_atoms(self, num_atoms_):
        self.num_atoms = num_atoms_
    
    def get_atom_types(self):
        return self.atom_types

    def set_num_atoms(self, atom_types_):
        self.atom_types = atom_types_
    
    def get_initial_positions(self):
        return self.initial_positions

    def set_initial_positions(self,initial_positions_):
        self.initial_positions = initial_positions_
    
    def get_positions(self):
        return self.positions

    def set_positions(self, positions_):
        self.positions = positions_
    
    def get_positions_store(self):
        return self.positions_store

    def set_positions_store(self, positions_store_):
        self.positions_store = positions_store_

    def assign_atom_types(self, input_file_):
        """
        Assign atom types.
        """
        data_file = get_data_file(input_file_)
        print "Obtaining atom types from " + data_file
        with open(data_file) as f:
            content = f.readlines()
        types_array=get_atom_types(input_file_)
        self.num_atom_types = int(str.split(content[3])[0])
        print str(self.num_atom_types) + " found of type:"
        print types_array
        for i in range(self.num_atoms):
            atom_type = int(str.split(content[14+(self.num_atom_types)+i])[1])
            self.atom_types.append(types_array[atom_type-1])
        print "Atom types are as follows:"
        print self.atom_types
        print "\n"


    
    def get_random_atom(self):
        """
        Selects which atom to move (drawn unifromally).
        """
        gen_rand = np.random.rand
        atom = int(self.num_atoms*gen_rand())
        return atom
    
    
    def move_atom_freely(self, atom_, radius_):
        """
        Moves atom_ to a randomly chosen point.
        The chosen point is within a  box of side r_.
        The box is centered at the original position x_[atom_] of atom_.
        Coords of atom_ are given by x_ (lammps object).
        """
        gen_rand = np.random.rand
        coord = gen_rand(3)-0.5
        self.positions[atom_][0] += radius_*coord[0]
        self.positions[atom_][1] += radius_*coord[1]
        self.positions[atom_][2] += radius_*coord[2]
    
    
    def move_atom_constrained(self, atom_, radius_):
        """
        Moves atom_ to a randomly chosen point.
        The chosen point is within a  box of side r_.
        The box is centered on the initial_position[atom_] of atom_.
        Coords of atom_ are given by x_ (lammps object).
        """
        if radius_ > 1.435: # Currently set to atomic sep
            radius_ = 1.435
        gen_rand = np.random.rand
        coord = gen_rand(3)-0.5
        self.positions[atom_][0] = self.initial_positions[3*atom_] + radius_*coord[0]
        self.positions[atom_][1] = self.initial_positions[3*atom_+1] + radius_*coord[1]
        self.positions[atom_][2] = self.initial_positions[3*atom_+2] + radius_*coord[2]


    def move_atom(self, atom_, radius_):
        """
        Dictates how atoms are moved.
        They can either move anywhere in simulation cell.
        Or be constrained within a given distance of their initial position.
        """
        if self.atom_types[atom_] == "Free":
            return
            self.move_atom_constrained(atom_, radius_)
        else:
            self.move_atom_freely(atom_, radius_) 
    

    def generate_configuration(self, lmp_, radius_,lower_bound_=float('-inf'),
                               upper_bound_=float('inf'), moves_ = 1000):
        """
        Generate an initial configuration that lies within an energy range.
        """
    
        # Function shortcuts for speed.
        move = self.move_atom
        pick_atom = self.get_random_atom
    
        # Initialise lammps for initial input.
        lmp_.command("run 0 post no")
        lmp_.command("variable e equal pe")
        # Extract the atomic coordinates.
        # This is an array of pointers.
        # Note: Changing these values will alter the values in lammps.
        self.positions = lmp_.extract_atom("x", 3)
    
        # Calculate and store the total energy of the initial configuration.
        lmp_.command("variable elast equal $e")
        elast = lmp_.extract_variable("elast", None, 0)
        
        # Store previous positions.
        self.positions_store = pos_to_array(lmp_) 
    
        for i in xrange(1, moves_):
            # Move all atoms
            # Gibbs sampler (all atoms are moved)
            for j in xrange(self.num_atoms):
                # Move selected atom to generate a new configuration.
                move(j, radius_)
            
            # Run lammps for the new configuration.
            lmp_.command("run 0 post no")
            # Refresh positions from lammps.
            # Required due to book-keeping (neighbour lists etc).
            self.positions = lmp_.extract_atom("x", 3)
            # Extract the total energy of the new configuration.
            lmp_.command("variable enew equal v_e")
            enew = lmp_.extract_variable("enew", None, 0)
            # Apply WL criterion to see if the move is accepted.
            if enew < upper_bound_ and enew > lower_bound_:
        #        print "Found initial configuration."
        #        print "Energy of initial configuration = " + str(enew)
                return 1./float(i)
            else:
                # If the move is rejected revert to original position of atom.
                array_to_pos(lmp_,self.positions_store)
        array_to_pos(lmp_,self.positions_store)
        lmp_.command("run 0 post no")
        #sys.exit("Failed to find configuration in energy range.")

