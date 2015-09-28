"""
Lammps MC class
"""
import os.path 
from lammps import lammps
from atom_commands import atoms

class lammps_system:
    """ A system defined within the Lammps MD program """


    def __init__(self, input_file_ = None, verbose_ = True):
        if input_file_ == None:
            print "Lammps input file not defined."
            exit()
        else:
            print "Using "+input_file_
        if verbose_:
            self.lmp = lammps()
        else:
            lammps_commands = ["-screen", "none", "-log", "none"]
            self.lmp = lammps("", lammps_commands)
        if os.path.isfile(input_file_):
            self.lmp.file(input_file_)
        else:
            print "Cannot find input file "+input_file_
        self.atom_system = atoms(self.lmp, input_file_ = input_file_)

        self.lmp.command("run 0 post no")
        self.lmp.command("variable e equal pe")
        self.atom_system.positions = self.lmp.extract_atom("x", 3)
    
    
    def get_energy(self):
        self.lmp.command("run 0 post no")
        self.atom_system.positions = self.lmp.extract_atom("x", 3)
        self.lmp.command("variable enew equal v_e")
        enew = self.lmp.extract_variable("enew", None, 0)
        return enew

    def generate_config(self):
        self.atom_system.generate_configuration(self.lmp, radius, moves_ = 1)

    def get_config(self):
        config = pos_to_array(self.lmp)

    def copy_config(self):


