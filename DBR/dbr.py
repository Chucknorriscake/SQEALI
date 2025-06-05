import sys
import meep as mp
import numpy as np
import os
import matplotlib.pyplot as plt 
import json
sys.path.append('../resonator_abstraction')
from meep_resonator import DBR

dbr_reference = DBR(config_file="dbr_reference.json")
pt = mp.Vector3(dbr_reference.base_vars["dbr_to_flux_2_distance"]+dbr_reference.base_vars["x_offset"],0)
dbr_reference.sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, pt, 1e-3))

reference_reflected_straight = dbr_reference.sim.get_flux_data(dbr_reference.fluxes[0])
straight_tran_flux = mp.get_fluxes(dbr_reference.fluxes[1])

dbr = DBR(config_file="dbr_main.json")
dbr.sim.load_minus_flux_data(dbr.fluxes[0], reference_reflected_straight)
dbr.sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, pt, 1e-3))

dbr_refl_flux = mp.get_fluxes(dbr.fluxes[0])
dbr_tran_flux = mp.get_fluxes(dbr.fluxes[1])

flux_freqs = mp.get_flux_freqs(dbr.fluxes[0])