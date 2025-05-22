# Linear Waveguide

- [Overview](#overview)
- [Configuration](#configuration)
    - [Base Variables](#base-variables)
    - [Derived Variables](#derived-variables)
    - [Geometry](#geometry)
    - [Source](#source)
    - [Flux](#flux)
- [Example Script](#example-script)
- [Documentation](#documentation)

## Overview
This example simulates coupling between two linear waveguides using Meep. The simulation analyzes the electric field distribution and calculates the coupling length based on the minima of the field envelope.

## Configuration

### Base Variables
```json
{
    "base_vars": {
        "wavelength": 3.0,
        "df": 0.2,
        "nfreq": 1000,
        "n": 1.8,
        "resolution": 16,
        "runtime": 500,
        "padding": 2,
        "dpml": 2,
        "waveguide_distance": 0.7,
        "waveguide_width": 0.5,
        "sx": 50
    },
    "derived_base_vars": {
        "fcen": "1/wavelength",
        "freqs_min": "fcen-df/2",
        "freqs_max": "fcen+df/2",
        "sy": "2*(waveguide_width + padding + dpml) + waveguide_distance"
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
    "coupling_fit_threshold": 0.15,
    "coupling_fit_first": 1,
    "coupling_fit_last": -1
}
```



### Geometry
```json
[
  {
    "type": "mp.Block",
    "size_x": "1e20",
    "size_y": "waveguide_width",
    "size_z": "1e20",
    "center_x": "0",
    "center_y": "waveguide_distance/2",
    "center_z": "0",
    "material_index": "n"
  },
  {
    "type": "mp.Block",
    "size_x": "1e20",
    "size_y": "waveguide_width",
    "size_z": "1e20",
    "center_x": "0",
    "center_y": "-waveguide_distance/2",
    "center_z": "0",
    "material_index": "n"
  }
]
```

### Source
```json
[
  {
    "usage": "simulation",
    "type": "mp.ContinuousSource",
    "component": "mp.Ez",
    "center_x": "-sx/2 + dpml+ padding/2",
    "center_y": "waveguide_distance/2",
    "center_z": "0",
    "size_x": "0",
    "size_y": "0.5*waveguide_width",
    "size_z": "0"
  }
]
```

### Flux
```json
[
  {
    "size_x": "0",
    "size_y": "0.5*waveguide_width",
    "size_z": "0",
    "center_x": "sx/2 - dpml - padding/2",
    "center_y": "waveguide_distance/2",
    "center_z": "0"
  },
  {
    "size_x": "0",
    "size_y": "0.5*waveguide_width",
    "size_z": "0",
    "center_x": "sx/2 - dpml - padding/2",
    "center_y": "-waveguide_distance/2",
    "center_z": "0"
  }
]
```


## Example script

```python
from meep_resonator import LinearWaveguides

# load and run the simulation
sim = LinearWaveguides("config_file.json")
sim.run_simulation(save_ez=True)

# plot the coupling lenghts
# by chaging threshold_percent=None and filter=[] change the fit parameters used

sim.plot_coupling_length(threshold_percent=0.1, filter=[0,-1])
# after the minima have been set, get the coupling length by running
coupling_length = sim.calc_coupling_length()
```

## Documentation
::: resonator_abstraction.meep_resonator.LinearWaveguides
