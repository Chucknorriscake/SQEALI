from meep_resonator import *
import numpy as np
import json
import timeit
import os



for conf_path in os.listdir(os.getcwd()):
    if  ("simulation_" in conf_path) and (conf_path != "simulation_sources.json") and (".json" in conf_path):
        ez_path = "ez_" + conf_path.split(".")[0] + "." + conf_path.split(".")[1] + ".csv"
        resonator = RingResonator(config_file=conf_path)
        resonator.run(resonator.run_resonance,ez_path, True) 
