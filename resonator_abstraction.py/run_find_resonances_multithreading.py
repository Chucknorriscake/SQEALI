from meep_resonator import *
import numpy as np

import timeit

start = timeit.default_timer()

wavelengths = np.linspace(1.5,2.5,16)
for wavelength in wavelengths:
    resonator = RingResonator(config_file="conf_overcoupling.json")
    resonator.base_vars["wavelength"] = wavelength
    resonator.run(resonator.find_resonances, resonance_conf_path="simulation_{}.json".format(wavelength)) 
    
    

stop = timeit.default_timer()
print('Time: ', stop - start)  
