from harmonic_oscillators import oscillator_system
from quantum_system import path_integral_system
from minimise import minimise
from wang_landau import wang_landau
number_os = 3
k_constant = 1.
hs = oscillator_system(number_ = number_os, k_constant_=k_constant)
test = path_integral_system(beads_ = 3, oscillator_ = hs)

wang_landau = wang_landau(test, dimensions_=2)
wang_landau.wl_monte_carlo(100)
#min_method.min_(10000000)
