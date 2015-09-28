#!/usr/bin/python
"""
Performs Wang-Landau Monte Carlo using Lammps as an energy solver.
"""

import os
import subprocess
import numpy as np
import copy

class wang_landau:
    """
    Wang Landau Class
    """
    def __init__(self, system_ = None, e_min_ = None, e_max_ = None,\
                e_min1_ = None, e_max1_ = None,\
                 range_ = None, lnf_ = None, dimensions_ = None ):

        if system_ == None:
            print "Please Define System"
        else:
            self.system = system_

        if dimensions_ == None:
            print "Using one-dimensional energy integral"
            self.dimension = 1
        else:
            self.dimension = dimensions_
            print "Using "+str(self.dimension)+"-dimensional energy integral"


        if lnf_ == None:
            print "Using default lnf"
            self.lnf = 1.
        else:
            self.lnf = lnf_
        
        self.wl_its = 500000000
        if self.dimension == 1: 
            if e_min_ == None and e_max_ == None:
                print "Using default energy range"
                self.e_min = 0.
                self.e_max = 50.
            elif e_min_ == None:
                print "Using default value for min energy"
                self.e_min = 0.
                self.e_max = e_max_
            elif e_max_ == None:
                print "Using default value for max energy"
                self.e_min = e_min_
                self.e_max = 50.
            else: 
                self.e_min = e_min_
                self.e_max = e_max_
                
            if range_ == None:
                print "Using default range"
                self.h_range = 20
            else:
                self.h_range = range_
            self.spacing = (self.e_max - self.e_min)/self.h_range
            self.histogram = np.zeros(self.h_range)
            self.dos = np.zeros(self.h_range)
        if self.dimension == 2:
            if e_min_ == None and e_max_ == None \
                    and e_min1_ == None and e_max1_ == None:
                print "Using default energy range"
                self.e_min = np.array([0.,0.])
                self.e_max = np.array([1.6,30.])
            elif e_min_ == None and e_max_ == None \
                    and e_min1_ == None and e_max1_ == None:
                print "I'm lazy please define whole range"
                exit()
            else: 
                self.e_min = np.array([e_min_,e_min1_])
                self.e_max = np.array([e_max_,e_max1_])
                
            if range_ == None:
                print "Using default number of histogram bars"
                self.h_range = 100
            else:
                self.h_range = range_
            self.histogram = np.zeros([self.h_range,self.h_range])
            self.dos = np.zeros([self.h_range,self.h_range])
            self.spacing = (self.e_max - self.e_min)/self.h_range

    def check_lnf(self):
        """
        Checks if each histogram region has been visited 1/sqrt(lnf) times.
        Used as a convergence criterion for WL.
        Note: This is the default convergence function.
        """
        min_val = np.amin(self.histogram[self.histogram > 0])
        if min_val > 1./np.sqrt(self.lnf):
            return True
        else:
            return False


    def check_flatness(self):
        """
        Checks if a histogram is flat to a given tolerance.
        Used as a convergence criterion for WL.
        This is the original Wang-Landau criterion.
        Note: This is not default must be uncommented to use.
        """
        mean = np.average(self.histogram[self.histogram > 0])
        min_val = np.amin(self.histogram[self.histogram > 0])
        if min_val > 0.8*mean: #0.8 is a magic number (i.e., completely arbitrary)
            return True
        else:
            return False

    def wl_monte_carlo(self, wl_total_):
        for i in range(wl_total_):
            result = self.wang_landau_simulation()
            if result:
                self.lnf = self.lnf/2.
                # this needs tidying up
                if self.dimension == 1:
                    self.histogram = np.zeros(self.h_range)
                elif self.dimension ==2:
                    self.histogram = np.zeros([self.h_range,self.h_range])
            else:
                print "Simulation failed to converge at iteration "+str(i)
                exit()

    def wang_landau_simulation(self):
        """
        Performs a Wang - Landau simulation.
        The WL simulation is carried out with update value lnf.
        """
        # Define commands for optimisation purposes.
        gen_rand = np.random.rand
        check_converge = self.check_lnf
        accept = 0 
        # Open path for Wang-Landau output.
        root_folder = "wl_output"
        if not os.path.exists(root_folder):
            os.mkdir(root_folder)
        path = root_folder+"/wl_output.txt"
        wl_output = open(path, 'a')
        
        if self.dimension == 1:
            elast = self.system.get_energy()      
            print "Energy of initial configuration = " + str(elast)
            wl_output.write("Energy of initial configuration = " + str(elast)+"\n")
            old_level = int((elast - self.e_min)/self.spacing)
        elif self.dimension == 2:
            elast = [self.system.get_energy1(),self.system.get_energy2()]
            old_level = (elast - self.e_min)/self.spacing
            old_level = [int(old_level[0]),int(old_level[1])]


    
        # Wang-Landau Monte Carlo Algorithm.
        for i in xrange(1, self.wl_its):
            # Store previous system.
            copy_system = copy.deepcopy(self.system.get_config())
            # Generate a new system 
            self.system.generate_config()
            enew = self.system.get_energy()
            # Extract the total energy of the new configuration.
            if self.dimension == 1:
                enew = self.system.get_energy()
                level = int((elast - self.e_min)/self.spacing)
                if enew < self.e_max and enew > self.e_min:
                    if(self.dos[level]-self.dos[old_level]) > 0.:
                        # If the move is accepted update stored energy.
                        elast = enew
                        old_level = level
                        # update acceptance counter
                        accept += 1

                    elif gen_rand() < \
                        np.exp(self.dos[level]-self.dos[old_level]):
                        # If the move is accepted update stored energy.
                        elast = enew
                        old_level = level
                        # update acceptance counter
                        accept += 1
                    else:
                        # If the move is rejected revert to original position of atom.
                        self.system.copy_config(copy_system) 
                    # Update histogram in energy space recording final configuration.
                    self.histogram[level] += 1
                    self.dos[level] += self.lnf
                else:
                    # If the move is rejected revert to original position of atom.
                    self.system.copy_config(copy_system) 
            elif self.dimension == 2:
                enew = [self.system.get_energy1(),self.system.get_energy2()]
                level = (elast - self.e_min)/self.spacing
                level = [int(level[0]),int(level[1])]
                if enew[0] < self.e_max[0] and enew[1] > self.e_min[0] and \
                enew[1] < self.e_max[1] and enew[1] > self.e_min[1]:
                    if gen_rand() < \
                        np.exp(self.dos[level[0],level[1]]-self.dos[old_level[0],old_level[1]]):
                        # If the move is accepted update stored energy.
                        elast = enew
                        old_level = level
                        # update acceptance counter
                        accept += 1
                    else:
                        # If the move is rejected revert to original position of atom.
                        self.system.copy_config(copy_system) 
                    # Update histogram in energy space recording final configuration.
                    self.histogram[level[0],level[1]] += 1
                    self.dos[level[0],level[1]] += self.lnf
                else:
                    # If the move is rejected revert to original position of atom.
                    self.system.copy_config(copy_system) 
            # Check convergence
            
            if (i > 1./np.sqrt(self.lnf)*self.h_range*self.h_range) and (i % 100 == 0) and i > 0:
                if check_converge():
                    print "Converged after " + str(i+1)+" iterations."
                    wl_output.write("Converged after " + str(i+1)+" iterations.\n")
                    print "Acceptance = " + str(float(accept)/float(i))+"\n"
                    wl_output.write("Acceptance = " \
                                    + str(float(accept)/float(i))+"\n")
                    path = root_folder+'/wl_dos_'+str(self.lnf)+'.txt'
                    if self.dimension==1:
                        e_range = np.zeros(self.h_range)
                        for i in range(self.h_range):
                            e_range[i] = self.e_min + i*self.spacing
                        write_dos = open(path, 'w')
                        for i in range(len(self.histogram)):
                            write_dos.write(str(e_range[i]) +
                                            '    ' + str(self.dos[i]) +
                                            '    ' + str(self.histogram[i]) +
                                            '\n')
                    elif self.dimension == 2:
                        write_dos = open(path, 'w')
                        e_range1 = np.zeros(self.h_range)
                        e_range2 = np.zeros(self.h_range)
                        for i in range(self.h_range):
                            e_range1[i] = self.e_min[0] + i*self.spacing[0]
                            e_range2[i] = self.e_min[1] + i*self.spacing[1]
                        for i in range(self.h_range):
                            for j in range(self.h_range):
                                write_dos.write(str(e_range1[i])+"   "+str(e_range2[j])+"   "+str(self.dos[i][j]) + "   \n")

                    return True

        return False
