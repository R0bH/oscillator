from harmonic_oscillators import oscillator_system
from quantum_system import path_integral_system
from minimise import minimise
number_os = 20
k_constant = 1.
hs = oscillator_system(number_ = number_os, k_constant_=k_constant)
test = path_integral_system(beads_ = 10, oscillator_ = hs)
for i in range(10):
    print test.ring_system[i].harmonic_oscillators
print  test.get_energy()
min_method = minimise(test)
min_method.min_(100000)
