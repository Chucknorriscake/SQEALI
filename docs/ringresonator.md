# Ring Resonator

- [Overview](#overview)
- [Configuration](#configuration)
    - [Stage 1](#stage-1)
        - [Base Variables](#base-variables)
        - [Harminv Source](#harminv-source)
    - [Stage 2](#stage-2)
        - [Generated Base Variables](#generated-base-variables)
        - [Continuous Source](#continuous-source)
    - [Stage 1 & 2](#stage-1-and-2)
        - [Geometry](#geometry)
        - [Flux](#flux)
- [Example Script](#example-script)
- [Documentation](#documentation)

## Overview
This simulation models a ring resonator coupled to a waveguide. It proceeds in two stages:
1. Estimate the resonance frequency using a Gaussian source and Harminv.
2. Simulate steady-state behavior using a continuous source at the resonance frequency.

## Configuration

The configuration is a twofold process. For step 1 we use a Gaussian source. and a wide range of frequencies to find resonances within the ring. For step two the script creates a new file that uses a continouus source and changes the wavelength to the resonance wavelength.

### Stage 1
#### Base Variables
```json
{
    "base_vars": {
        "wavelength": 3.0,
        "df": 0.2,
        "nfreq": 1000,
        "radius": 2,
        "ring_width": 0.4,
        "n": 1.8,
        "resolution": 16,
        "runtime": 5000,
        "padding": 2,
        "dpml": 2,
        "waveguide_padding": 0.23
    },
    "derived_base_vars": {
        "fcen": "1/wavelength",
        "freqs_min": "fcen-df/2",
        "freqs_max": "fcen+df/2",
        "outer_radius": "radius + ring_width",
        "waveguide_width": "ring_width",
        "sy": "2*(outer_radius + padding + dpml)",
        "waveguide_xpos": "-(radius+ ring_width + waveguide_padding)",
        "sx": "sy + 50"
    },
    "cell": {
        "size_x": "sx",
        "size_y": "sy",
        "size_z": "0"
    },
    "boundary_layer": [
        {
            "type": "mp.PML",
            "depth": "dpml"
        }
    ],
    "geometry_file": "geometry.json",
    "flux_file": "flux.json",
    "source_file": "harminv_sources.json",

}
```





#### Harminv Source
```json
[
  {
    "usage": "harminv",
    "type": "mp.GaussianSource",
    "component": "mp.Ez",
    "center_x": "-sx/2 + dpml + 0.5",
    "center_y": "waveguide_xpos",
    "center_z": "0",
    "size_x": "0",
    "size_y": "0.5*ring_width",
    "size_z": "0"
  }
]
```

### Stage 2 
These file is generated automatically.
#### Generated Base Variables
```json
{
    "base_vars": {
        "wavelength": "[resonance from step 1]",
        "df": 0.02,
        "nfreq": 1000,
        "radius": 2,
        "ring_width": 0.4,
        "n": 1.8,
        "resolution": 16,
        "runtime": 2500,
        "padding": 2,
        "dpml": 2,
        "waveguide_padding": 0.23
    },
    "derived_base_vars": {
        "fcen": "1/[resonance from step 1]",
        "freqs_min": "fcen-df/2",
        "freqs_max": "fcen+df/2",
        "outer_radius": "radius + ring_width",
        "waveguide_width": "ring_width",
        "sy": "2*(outer_radius + padding + dpml)",
        "waveguide_xpos": "-(radius+ ring_width + waveguide_padding)",
        "sx": "sy + 50"
    },
    "cell": {
        "size_x": "sx",
        "size_y": "sy",
        "size_z": "0"
    },
    "boundary_layer": [
        {
            "type": "mp.PML",
            "depth": "dpml"
        }
    ],
    "source_file": "simulation_sources.json",
    "geometry_file": "geometry.json",
    "flux_file": "flux.json",
    "ez_file": "ez_simulation_[resonance from step 1].csv"
}
```
#### Continuous Source
```json
[
  {
    "usage": "simulation",
    "type": "mp.ContinuousSource",
    "component": "mp.Ez",
    "center_x": "-sx/2 + dpml + 0.5",
    "center_y": "waveguide_xpos",
    "center_z": "0",
    "size_x": "0",
    "size_y": "0.5*ring_width",
    "size_z": "0"
  }
]
```


### Stage 1 and 2

Both stages use the same geometry and flux configuration.
#### Geometry
```json
[
    {
        "type": "mp.Cylinder",
        "radius": "outer_radius",
        "center_x": "0",
        "center_y": "0",
        "center_z": "0",
        "material_index": "n"
    },
    {
        "type": "mp.Cylinder",
        "radius": "radius",
        "center_x": "0",
        "center_y": "0",
        "center_z": "0",
        "material_index": "1"
    },
    {
        "type": "mp.Block",
        "size_x": "1e20",
        "size_y": "waveguide_width",
        "size_z": "1e20",
        "center_x": "0",
        "center_y": "waveguide_xpos",
        "center_z": "0",
        "material_index": "n"
    }
]
```

#### Flux
```json
[
    {
        "size_x": "0",
        "size_y": " 0.5*ring_width",
        "size_z": "0",
        "center_x": "-sx/2 +dpml + 0.5+3",
        "center_y": "waveguide_xpos",
        "center_z": "0"
    },
    {
        "size_x": "0",
        "size_y": " 0.5*ring_width",
        "size_z": "0",
        "center_x": "sx/2 -dpml - 0.5",
        "center_y": "waveguide_xpos",
        "center_z": "0"
    }
]
```

## Example Script

This two-step simulation process finds resonant frequencies and then runs full simulations for each using Meep and the `RingResonator` class.

---

### Step 1: Sweep Over Wavelengths and Find Resonances

This script updates the configuration file for each wavelength, initializes a `RingResonator` object, and finds resonant frequencies using `find_resonances`.

```python
from meep_resonator import *
import numpy as np
import json
import timeit

start = timeit.default_timer()

# Define wavelength sweep range
wavelengths = np.linspace(0.6, 10, 30)

# Path to base configuration file
conf_path = "conf_overcoupling.json"

for wavelength in wavelengths:
    try:
        # Load base config and update wavelength
        with open(conf_path) as config_file:
            config = json.load(config_file)
            config["base_vars"]["wavelength"] = wavelength

        # Save updated config
        with open(conf_path, 'w') as f:
            json.dump(config, f, indent=4)

        # Initialize resonator and run harminv-based frequency detection
        resonator = RingResonator(config_file=conf_path)
        print(f"Running resonance detection at wavelength: {resonator.base_vars['wavelength']}")
        resonator.find_resonances(
            resonance_conf_path=f"simulation_{wavelength}.json"
        )
        
    except Exception as e:
        raise FileNotFoundError(f"Error processing wavelength {wavelength}: {e}")

stop = timeit.default_timer()
print('Total Time: ', stop - start)
```

---

### Step 2: Run Simulations for Each Resonance Config

This script loops through all generated resonance config files and runs the actual simulation using the detected resonant frequency.

```python
from meep_resonator import *
import os

# Iterate over all simulation config files
for conf_path in os.listdir(os.getcwd()):
    if ("simulation_" in conf_path) and (conf_path != "simulation_sources.json") and (conf_path.endswith(".json")):
        print(f"Running simulation with config: {conf_path}")
        resonator = RingResonator(config_file=conf_path)
        resonator.run(resonator.run_resonance, store_ez=True)
```

## Documentation

::: resonator_abstraction.meep_resonator.RingResonator