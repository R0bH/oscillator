from harmonic_oscillators import oscillator_system
from quantum_system import path_integral_system
from lammps_system import lammps_system
from minimise import minimise
from wang_landau import wang_landau
number_os = 3
k_constant = 1.
#hs = oscillator_system(number_ = number_os, k_constant_=k_constant)
lp = lammps_system(input_file_ = "iron.in")
test = path_integral_system(beads_ = 3, system_ = lp)
print lp.get_energy()
print test.get_energy()
wang_landau = wang_landau(test, dimensions_=2,e_min_ = -10., e_max_=10., e_min1_ = 0., e_max1_=100.,range_=100)
wang_landau.wl_monte_carlo(100)
#min_method.min_(10000000)
