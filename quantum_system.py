import numpy as np
import copy 

from harmonic_oscillators import oscillator_system

class path_integral_system:
    """
    N Bead Path Integral Formalism
    Defines N systems coupled by harmonic beads
    Mehods are named to allow polymorphism.
    """
    def __init__(self, beads_ = None, oscillator_ = None):
        if beads_ == None:
            print "No beads defined."
            exit()
        elif oscillator_ == None:
            print "No oscillator system to clone."
            exit()
        else:
            self.beads = beads_
            self.ring_system = []
            for i in range(beads_):
                self.ring_system.append(oscillator_system(number_ = len(oscillator_.get_config()), k_constant_ = oscillator_.get_k_constant()))

    def get_config(self):
        return self.ring_system

    def set_config(self, ring_system_):
        self.ring_system = ring_system_
    
    def copy_config(self, ring_system_):
        self.ring_system = copy.deepcopy(ring_system_)

    def generate_config(self):
        for i in range(self.beads):
            self.ring_system[i].generate_config()

    def get_energy(self):
        """
        Calculates total energy of the system.
        Potential energy is average of the sub systems.
        Oscillator energy: each oscillator in a classic system interacts
        with its counterpart in i-1 and i+1 system.
        """
        # Calculate potential of systems
        pot = 0.
        for i in range(self.beads):
            pot += self.ring_system[i].get_energy()

        # calculate energy of ring
        osc = 0.
        for i in range(self.beads):
            if i == self.beads-1:
                for j in range( len(self.ring_system[0].get_config())):
                    osc += (self.ring_system[i].get_config_element(j) - self.ring_system[0].get_config_element(j))**2
            else: 
                for j in range( len(self.ring_system[0].get_config())):
                    osc += (self.ring_system[i].get_config_element(j) - self.ring_system[i+1].get_config_element(j))**2
        return pot/self.beads + osc

    def get_potential_energy(self):
        """
        Calculates potential energy of the system.
        Potential energy is average of the sub systems.
        """
        # Calculate potential of systems
        pot = 0.
        for i in range(self.beads):
            pot += self.ring_system[i].get_energy()

        return pot/self.beads

    def get_oscillator_energy(self):
        """
        Calculates oscillator energy of the system.
        Oscillator energy: each oscillator in a classic system interacts
        with its counterpart in i-1 and i+1 system.
        """
        # calculate energy of ring
        osc = 0.
        for i in range(self.beads):
            if i == self.beads-1:
                for j in range( len(self.ring_system[0].get_config())):
                    osc += (self.ring_system[i].get_config_element(j) - self.ring_system[0].get_config_element(j))**2
            else: 
                for j in range( len(self.ring_system[0].get_config())):
                    osc += (self.ring_system[i].get_config_element(j) - self.ring_system[i+1].get_config_element(j))**2
        return osc

