
# Advanced Photonics MEEP Simulation

## Description
This framework streamlines the creation of MEEP simulations, transitioning from a code-based approach to an easily editable configuration file approach, where all simulation parameters are handled through `.json` files.

## Installation
Prerequisites for this framework include a running version of Anaconda or Miniconda installed on your computer.

To install the environment, run:

```bash
conda create --name <env> --file requirements.txt
conda activate <env>
```

The framework was developed and tested in `Python 3.11.11`. The main focal packages are:

- pymeep (multithreaded)
- numpy
- scipy
- matplotlib

## Ring resonator: quick guide
For single threading:
``` python
from meep_resonator import *

# run find resonance

resonator = RingResonator("find_resonance.json")
resonance.find_resonances("resonance.json")

# run simulation on resonance

resonator = RingResonator("resonance.json")
resonantor.run_resonances

```
For multithreading:
```python 
from meep_resonator import *

# run find resonance

wavelengths = np.linspace(0.6,10,30)
conf_path = "find_resonance.json"

for wavelength in wavelengths:
    try:
        with open(conf_path) as config_file:

            config = json.load(config_file)
            config["base_vars"]["wavelength"] = wavelength

        with open(conf_path, 'w') as f:
            json.dump(config, f, indent=4)
        resonator = RingResonator(config_file=conf_path)
        resonator.run(resonator.find_resonances, resonance_conf_path="simulation_{}.json".format(wavelength)) 
    except:
        raise FileNotFoundError 



# run simulation on resonance

for conf_path in os.listdir(os.getcwd()):
    if  ("simulation_" in conf_path) and (conf_path != "simulation_sources.json") and (".json" in conf_path):
        ez_path = "ez_" + conf_path.split(".")[0] + "." + conf_path.split(".")[1] + ".csv"
        resonator = RingResonator(config_file=conf_path)
        resonator.run(resonator.run_resonance,ez_path, True) 


```
## parallel waveguides: quick guide
### DISCLAIMER: not implimented

## Usage

Using this library is straightforward. Once the configuration files, described in **Config Files: Simulation Prerequisites**, are set up, the resonator can be created. The configuration file path is passed to the simulation as a string. Currently, the simulation options offer two paths inheriting their setup properties from the `SimulationBase` class.

### Ring Resonators
To use ring resonators:

```python
resonator = RingResonator("[path to main config.json]")
```

This loads the desired ring resonator configuration into the class and automatically sets up the MEEP environment.  
Since the resonance frequency of this setup is originally unknown, follow these steps:

1. Perform a Harminv simulation with a Gaussian source by running:
   ```python
   resonator.find_resonances("path_to_result_config")
   ```
   This automatically creates a preconfigured, updated configuration file for step two. By default, the path to the source file is changed to `"simulation_source.json"`, as the provided simulation source switches to an `mp.ContinuousSource`.

2. Conduct a full-length simulation of the dynamics at the resonance frequencies using the new configuration files. Initialize a new resonator object with:
   ```python
   resonator = RingResonator("path_to_result_config")
   ```
   Then, run a normal simulation with:
   ```python
   resonator.run_resonances(save_ez=True)
   ```
   The parameter `save_ez` is used to save the `Ez` field data of a simulation to a CSV file, if needed.

### Parallel Waveguides

#### DISCLAIMER: Not implemented
```python
parallel = ParallelWaveguides("[path to main config.json]")
```

## Base Simulation
```python
simpleBase = SimulationBase("[path to main config.json]")
```
The `SimulationBase` class serves as the foundation, reading all the required configuration files and setting up the basic simulation environment.

## Config Files: Simulation Prerequisites

To run a simulation, this framework requires three main components:

- A main configuration file
- A source configuration file
- A geometry configuration file

An optional fourth file, a flux configuration, can be created to add fluxes to the simulation. However, it is not strictly necessary.

### DISCLAIMER:  
**Known Issue:** `mp.inf` currently does not work within the configuration files. However, MEEP sets `mp.inf` to `1e20`. As a workaround, substitute every `mp.inf` with `1e20`.

### Main Configuration

The main configuration JSON file manages all basic simulation variables. While most variable keys are fully customizable, three keys are mandatory within `base_vars` and `derived_base_vars` for the script to work: `wavelength`, `df`, and `fcen`.

For all other structures, the keys must remain as defined.

```json
{
    "base_vars": {          
        "wavelength": 1,
        "df": 0.1,
        "sx": 1
    },
    "derived_base_vars": {
        "fcen": "1/wavelength",
        "sy": "sx + 2"
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
    "source_file": "harminv_sources.json",
    "geometry_file": "geometry.json",
    "flux_file": "flux.json"
}
```

### Sources Configuration

The source configuration is a list of all sources to be used within the simulation domain. The keys align with MEEP documentation for source creation parameters. An additional key, `usage`, identifies the purpose of the source file.

```json
[
    {
        "usage": "harminv",
        "type": "mp.GaussianSource",
        "component": "mp.Ez",
        "center_x": "sx-sx/2",
        "center_y": "0",
        "center_z": "0",
        "size_x": "0",
        "size_y": "0.5",
        "size_z": "0"
    }
]
```

### Geometry Configuration

The geometry configuration is a list of geometries used in the simulation. Currently supported geometry types include `mp.Cylinder` and `mp.Block`. These configurations align with the MEEP documentation.

```json
[
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
        "size_y": "",
        "size_z": "1e20",
        "center_x": "0",
        "center_y": "waveguide_xpos",
        "center_z": "0",
        "material_index": "n"
    }
]
```

### Flux Configuration

If a flux configuration is needed for the simulation, the configuration file must follow this structure:

```json
[
    {
        "size_x": "0",
        "size_y": "0.5*ring_width",
        "size_z": "0",
        "center_x": "-sx/2 +padding/2 + 0.5+3",
        "center_y": "waveguide_xpos",
        "center_z": "0"
    },
    {
        "size_x": "0",
        "size_y": "0.5*ring_width",
        "size_z": "0",
        "center_x": "sx/2 -padding/2 - 0.5",
        "center_y": "waveguide_xpos",
        "center_z": "0"
    }
]
```

## Roadmap

## Author
Julian Verhey

## Project Status
This is a very crude work-in-progress repository. The codebase is, at best, spaghetti and on fire. The proposed workflow has not been thoroughly tested, and edge cases and error handling are not fully implemented. However, as this project is a fun proof of concept and an ode to OOP—as well as sinking more time into a wonderful programming task than necessary—I am perfectly okay with this and happy about it.

