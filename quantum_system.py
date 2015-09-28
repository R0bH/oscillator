import numpy as np
import copy 


class path_integral_system:
    """
    N Bead Path Integral Formalism
    Defines N systems coupled by harmonic beads
    Mehods are named to allow polymorphism.
    """
    def __init__(self, beads_ = None, system_ = None):
        if beads_ == None:
            print "No beads defined."
            exit()
        elif system_ == None:
            print "No system to clone."
            exit()
        else:
            self.beads = beads_
            self.ring_system = []
            for i in range(beads_):
                #self.ring_system.append(oscillator_system(number_ = len(oscillator_.get_config()), k_constant_ = oscillator_.get_k_constant()))
                system_.generate_config()
                clone = copy.deepcopy(system_)
                self.ring_system.append(clone)

    def get_config(self):
        return self.ring_system

    def set_config(self, ring_system_):
        self.ring_system = ring_system_
    
    def copy_config(self, ring_system_):
        self.ring_system = copy.deepcopy(ring_system_)

    def generate_config(self):
        bead_move = int(np.random.random()*self.beads)
#        for i in range(self.beads):
#            self.ring_system[i].generate_config()
        self.ring_system[bead_move].generate_config()

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
        osc1 = 0.
        for i in range(self.beads):
            if i == self.beads-1:
                osc += sum((self.ring_system[i].get_config() - self.ring_system[0].get_config())**2)
            else: 
                osc += sum((self.ring_system[i].get_config() - self.ring_system[i+1].get_config())**2)
        return pot/self.beads + 0.5*osc*self.beads

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
                osc += sum((self.ring_system[i].get_config() - self.ring_system[0].get_config())**2)
            else: 
                osc += sum((self.ring_system[i].get_config() - self.ring_system[i+1].get_config())**2)
        return 0.5*osc*self.beads
    
    # Alias methods for PIWL
    get_energy1 = get_potential_energy
    get_energy2 = get_oscillator_energy
