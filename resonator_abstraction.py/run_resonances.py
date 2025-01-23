from meep_resonator import *
import numpy as np
import json
import timeit
import os


start = timeit.default_timer()

for conf_path in os.listdir(os.getcwd()):
    if ("simulation_" in conf_path) and (conf_path != "simulation_sources.json"):
        ez_path = "ez_" + conf_path.split(".")[0] + "." + conf_path.split(".")[1] + ".csv"
        print(ez_path)
        resonator = RingResonator(config_file=conf_path)
        resonator.run(resonator.run_resonance,ez_path, True) 


stop = timeit.default_timer()
print('Time: ', stop - start)  