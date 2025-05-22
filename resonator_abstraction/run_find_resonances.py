from meep_resonator import *
import numpy as np

import timeit
import json
start = timeit.default_timer()

wavelengths = np.linspace(0.6,10,30)

conf_path = "conf_overcoupling.json"

for wavelength in wavelengths:
    try:
        with open("conf_overcoupling.json") as config_file:

            config = json.load(config_file)
            config["base_vars"]["wavelength"] = wavelength

        with open(conf_path, 'w') as f:
            json.dump(config, f, indent=4)
        
        resonator = RingResonator(config_file=conf_path)
        print(resonator.base_vars["wavelength"])
        resonator.base_vars["wavelength"] = wavelength
        resonator.find_resonances(resonance_conf_path="simulation_{}.json".format(wavelength))
    except:
        raise FileNotFoundError 
    
    
    

stop = timeit.default_timer()
print('Time: ', stop - start)  
