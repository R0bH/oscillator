import numpy as np
import copy

from harmonic_oscillators import oscillator_system

class minimise:
    """
    Minimise method class
    """
    def __init__(self, system_ = None):
        if system_ == None:
            print "Please define operator system"
            exit()
        else:
            self.system = system_
            self.min_system = copy.deepcopy(system_)

    def gen_config(self):
        self.system.generate_config()
        store = self.system.get_config()
        #print store[0].get_config()
        #print self.system.get_energy()
    def check_min(self):
        if self.system.get_energy() < self.min_system.get_energy():
            self.min_system.copy_config(self.system.get_config())
  
    def min_(self, its_):
        for i in range(its_):
            self.gen_config()
            self.check_min()
         #   print self.system.get_energy()
            print self.min_system.get_energy()
