"""
Commands used to select and move atoms.
"""

import numpy as np
import copy

class oscillator_system:
    """ A system of harmonic oscilators """
    def __init__(self, number_ = None,\
            k_constant_ = None, x_range_ = None):
        if x_range_ == None:
            print "Range not defined assumin +/- 1."
            self.x_range = 2.
        else:
            self.x_range = x_range_

        if (number_ == None) and (k_constant_ == None):
            self.harmonic_oscillators = []
            self.k_constant = []
        
        if (number_ == None):
            if (type(k_constant_) is int) \
                or (type(k_constant_) is float):
                    print "Please define number of oscillators"
                    exit()
            if (type(k_constant_) is list):
                self.harmonic_oscillators = self.x_range*np.random.random(len(k_constant_)) - self.x_range/2.
                self.k_constant = k_constant_
        
        if (k_constant_ == None):
            self.harmonic_oscillators = self.x_range*np.random.random(number_) - self.x_range/2.
            self.k_constant = k_constant_
        
        else:
            if (type(k_constant_) is int) \
                or (type(k_constant_) is float):
                self.harmonic_oscillators = self.x_range*np.random.random(number_) - self.x_range/2.
                self.k_constant = k_constant_
            elif number_ == len(k_constant_):
                self.harmonic_oscillators = self.x_range*np.random.random(number_) - self.x_range/2.
                self.k_constant = k_constant_
            else:
                print "Number of osciallators and number of constants do not match"
                exit()

    def get_config(self):
        return self.harmonic_oscillators

    def set_config(self, harmonic_oscillators_):
        self.harmonic_oscillators = harmonic_oscillators_
    
    def copy_config(self, harmonic_oscillators_):
        self.harmonic_oscillators = copy.deepcopy(harmonic_oscillators_)
    
    def get_config_element(self,element_):
        return self.harmonic_oscillators[element_]
    
    def set_config_element(self,element_, value_):
        self.harmonic_oscillators[element_] = value_
    
    def copy_config_element(self, harmonic_oscillators_):
        self.harmonic_oscillators = copy.deepcopy(harmonic_oscillators_)
    
    def get_k_constant(self):
        return self.k_constant

    def set_k_constant(self, k_constant_):
        self.k_constant = k_constant_
    
    def copy_k_constant(self, k_constant_):
        self.k_constant = copy.deepcopy(k_constant_)

    def get_x_range(self):
        return self.x_range

    def set_x_range(self, x_range_):
        self.x_range = x_range_

    def generate_config(self, oscillator_ = None):
        if oscillator_ == None:
            self.generate_harmonic_oscillators()
        elif oscillator_ == True:
            self.change_random_oscillator()
        else:
            print "oscillator_ improperly defined"
            exit()

    def generate_harmonic_oscillators(self, number_ = None):
        if number_ == None:
            self.harmonic_oscillators = \
                    self.x_range*np.random.random(len(self.harmonic_oscillators)) \
                    - self.x_range/2.
        else:
            self.harmonic_oscillators = self.x_range*np.random(number_) - self.x_range/2.
    
    def change_random_oscillator(self):
        oscillator = int(np.random.random()\
                *len(self.get_harmonic_oscillators()))
        self.harmonic_oscillators[oscillator] = \
                self.x_range*np.random.random() - self.x_range/2.

    def get_energy(self):
        """
        Energy of system given by 0.5*k**2
        May make oscillators coupled but not sure.
        May also make 2d or 3d
        """
        if self.k_constant == None:
            print "k_constant not defined"
            exit()
        elif (type(self.k_constant) is int) \
                or (type(self.k_constant) is float):
            return 0.5*self.k_constant*np.sum(self.harmonic_oscillators**2)
        elif (type(self.k_constant) is list) and (len(self.harmonic_oscillators) == len(self.k_constant)):
            return 0.5*np.sum(self.k_constant*self.harmonic_oscillators**2)
        else:
            print "k_constant improperly defined"
            exit()


