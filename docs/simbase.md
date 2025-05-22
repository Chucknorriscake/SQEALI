# Simulation Base

- [Overview](#overview)
- [Configuration](#configuration)
    - [Base](#base)
    - [Source](#source)
    - [Geometry](#geometry)
    - [Flux](#flux)
- [Documentation](#documentation)
## Overview

This class is not meant to be run standalone. It serves as a base creation layer for inherting child classes for simulations.
the Simulation base provides functionality in creating a MEEP simulation domain and creating geometries, sources, fluxes, the cell, boundary layers aswell as abstracting exact positioning to a simple base_vars dictionary.

## Configuration

### Base

All inherting classes expect (at minimum) the following structure from a configuration json file

| Key               | Description                         |
|------------------|-------------------------------------|
| `base_vars`       | Number-based base variables         |
| `derived_base_vars` | Calculated (derived) base variables |
| `cell`            | Size of the simulation cell         |
| `boundary_layer`  | Type and depth of the boundary layer |
| `source_file`  | Path to the sources |
| `geometry_file`  | Path to the geometry |
| `flux_file`  | Path to the flux |

```json
{
    "base_vars": {
        "base_var_1": "...",
        "base_var_2": "..."
    },
    "derived_base_vars": {
        "derived_base_var_1": "...",
        "derived_base_var_2": "..."
        
    },
    "cell": {
        "size_x": "...",
        "size_y": "...",
        "size_z": "..."
    },
    "boundary_layer": [
        {
            "type": "...",
            "depth": "..."
        }
    ],
    "source_file": "harminv_sources.json",
    "geometry_file": "geometry.json",
    "flux_file": "flux.json"
}
```

### Source

The sources file needs at least the type, component, location and size of a source. The type describes if it is a countinuous source or a gaussian source and follows the naming scheme of MEEP (mp.GaussianSource or mp.ContinouusSource). The component also follows MEEP (mp.Ez for source in z direction). Note that this json is a list. This is necessairy as multiple sources are allowed.

```json
[
    {   
        "type": "...",
        "component": "...",
        "center_x": "...",
        "center_y": "...",
        "center_z": "...",
        "size_x": "...",
        "size_y": "...",
        "size_z": "..."
    }
]
```

### Geometry

currently implemented are mp.Cylinders and mp.Blocks. Sizes and locations can be calculated by the usage of `base_vars` and `derived_base_vars`. Currently the use of `mp.inf` is not supported. Hence it is necessairy to use 1e20 (which is equivalent to mp.inf according to MEEP documentation) 
```json
[
    {
        "type": "mp.Cylinder",
        "radius": "...",
        "center_x": "...",
        "center_y": "...",
        "center_z": "...",
        "material_index": "..."
    },
    {
        "type": "mp.Block",
        "size_x": "1e20",
        "size_y": "...",
        "size_z": "1e20",
        "center_x": "...",
        "center_y": "...",
        "center_z": "...",
        "material_index": "..."
    }
]

```
### Flux

Fluxes need the following parameters
```json
[
    {
        "size_x": "...",
        "size_y": "...",
        "size_z": "...",
        "center_x": "...",
        "center_y": "...",
        "center_z": "..."
    }
]
```
## Documentation

::: resonator_abstraction.meep_resonator.SimulationBase