import numpy as np
import matplotlib.pyplot as plt
import meep as mp
import json
import logging
import threading
import queue
from IPython.display import Video
from scipy.optimize import curve_fit
from helpers import parse_expression
from tqdm import tqdm
import csv 
# Configure logging for debugging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Base class for handling configuration and setup
class SimulationBase:
    '''
    Meep wrapper that allows to build simulations by providing a framework that loads a JSON config file and sets up
    the simulation domain accordingly to the file. 
    '''
    def __init__(self, config_file : str) -> None:
        '''
        just an init
        
        Keywords:
            config_file: path of the config file that contains the simulation parameters
        '''
        self.base_vars = {}
        self.sources = []
        self.geometry = []
        self.boundary_layers = []
        self.simulation_domain = mp.Vector3(10, 10, 0)
        self.resolution = 16
        self.runtime = 3000
        self.config_file=config_file
        self.sim = None
        self.ez_field = None
        # parallel processing 
        mp.divide_parallel_processes(1)
        # shut up if needed 
        mp.verbosity.meep = 0
        # multithreading support
        self._worker_thread = None
        self._result_queue = queue.Queue()
        self._stop_event = threading.Event()
        # progress bar 
        self.progress_bar = None # lifetime is set directly before starting the simulation
        # load config file to setup the simulation
        if config_file:
            self.sim, self.sources_path ,self.flux_path,self.geometry_path = self._load_conf(config_file)
            self._build_flux(self.flux_path)
     
#################################################################
# Internal functions 
#################################################################

    def _load_conf(self, path: str) -> mp.Simulation:
        '''
        loads the configuration file from a path containing JSON and returns a setup Meep Simulation object
        
        Keywords:
            path: location of config file (JSON)
        '''
        try:  
            with open(path) as config_file:
                config = json.load(config_file)

                sources_path = config["source_file"]
                flux_path = config["flux_file"]
                geometry_path = config["geometry_file"]
                # load ez if file available:
                try:
                    ez_path = config["ez_file"]
                    self.ez_field = self._load_ez_data(ez_path)
                    logger.info("Simulation parameters have been run before and ez field has been loaded!")
                except:
                    pass
                    
            self.base_vars = self._load_basic_variables(path)
            self.sources = self._build_sources(sources_path)
            self.geometry = self._build_geometry(geometry_path)
            self.boundary_layers = self._build_boundary_layers(path)
            self.simulation_domain = self._build_cell(path)
            self._set_simulation_params()
            

            return mp.Simulation(
                cell_size=self.simulation_domain,
                boundary_layers=self.boundary_layers,
                geometry=self.geometry,
                sources=self.sources,
                resolution=self.resolution,
            ),sources_path, flux_path, geometry_path 
        except Exception as e:
            logger.error(f"Error loading configuration: {e}")
            raise

    def _load_basic_variables(self, path: str) -> dict:
        '''
        basic variables define the geometric parameters of a Meep simulation, as well as the sources frequency and width.
        these are loaded from the main config file directly
        
        Keywords:
            path: location of main config file
            
        Return:
            base_vars: dictionary of all basic variables from which the rest of the simulation domain will be set up 
        '''
        try:
            with open(path) as config_file:
                config = json.load(config_file)
                base_vars = config.get("base_vars", {})
                derived_vars = config.get("derived_base_vars", {})
                for key, val in derived_vars.items():
                    base_vars[key] = parse_expression(base_vars, val)

                base_vars["freqs"] = np.linspace(
                    base_vars["freqs_min"], base_vars["freqs_max"], base_vars["nfreq"]
                )
                return base_vars
        except KeyError as e:
            logger.error(f"Missing key in configuration: {e}")
            raise
        except Exception as e:
            logger.error(f"Error loading basic variables: {e}")
            raise

    def _set_simulation_params(self) -> None:
        '''
        Sets the resolution and runtime of the simulation domain
        '''
        self.resolution = self.base_vars.get("resolution", 16)
        self.runtime = self.base_vars.get("runtime", 100)

    def _build_sources(self, path: str) -> list:
        '''
        builds the sources that are defined within the sources file
        
        Keywords: 
            path: location of the sources file (JSON)
        Return:
            sources: list of mp.Sources that can be fed into mp.Simulation
        '''
        try:
            with open(path) as f:
                sources_config = json.load(f)
                sources = []

                for source_conf in sources_config:
                    source_type = eval(source_conf["type"])
                    component = eval(source_conf["component"])
                    center = mp.Vector3(
                        parse_expression(self.base_vars, source_conf["center_x"]),
                        parse_expression(self.base_vars, source_conf["center_y"]),
                        parse_expression(self.base_vars, source_conf["center_z"]),
                    )
                    size = mp.Vector3(
                        parse_expression(self.base_vars, source_conf["size_x"]),
                        parse_expression(self.base_vars, source_conf["size_y"]),
                        parse_expression(self.base_vars, source_conf["size_z"]),
                    )

                    sources.append(
                        mp.Source(
                            src=source_type(
                                frequency=self.base_vars["fcen"],
                                fwidth=self.base_vars["df"],
                            ),
                            component=component,
                            center=center,
                            size=size,
                        )
                    )
                return sources
        except Exception as e:
            logger.error(f"Error building sources: {e}")
            raise

    def _build_geometry(self, path: str) -> list:
        '''
        builds the geometries that are defined within the geometry file. Right now ONLY cylinders and blocks are supported.
        
        Keywords: 
            path: location of the geometries file (JSON)
        
        Return:
            geometries: list of mp.Cylinder, mp.Blocks that can be fed into mp.Simulation
        '''
        try:
            with open(path) as f:
                geometry_config = json.load(f)
                geometry = []

                for geom_conf in geometry_config:
                    center = mp.Vector3(
                        parse_expression(self.base_vars, geom_conf["center_x"]),
                        parse_expression(self.base_vars, geom_conf["center_y"]),
                        parse_expression(self.base_vars, geom_conf["center_z"]),
                    )
                    material_index = parse_expression(self.base_vars, geom_conf["material_index"])

                    if geom_conf["type"] == "mp.Block":
                        size = mp.Vector3(
                            parse_expression(self.base_vars, geom_conf["size_x"]),
                            parse_expression(self.base_vars, geom_conf["size_y"]),
                            parse_expression(self.base_vars, geom_conf["size_z"]),
                        )
                        geometry.append(
                            mp.Block(
                                size=size,
                                center=center,
                                material=mp.Medium(index=material_index),
                            )
                        )
                    elif geom_conf["type"] == "mp.Cylinder":
                        radius = parse_expression(self.base_vars, geom_conf["radius"])
                        geometry.append(
                            mp.Cylinder(
                                radius=radius,
                                material=mp.Medium(index=material_index),
                                center=center,
                            )
                        )
                return geometry
        except Exception as e:
            logger.error(f"Error building geometry: {e}")
            raise

    def _build_boundary_layers(self, path: str) -> list:
        '''
        creates the simulations boundary layer. Read directly from the config file
        
        Keywords: 
        path: location of the main config file
        
        Return:
            boundary_layers: list of boundary layers
        '''
        try:
            with open(path) as config_file:
                config = json.load(config_file)
                boundary_layers = []

                for layer_conf in config.get("boundary_layer", []):
                    layer_type = eval(layer_conf["type"])
                    depth = parse_expression(self.base_vars, layer_conf["depth"])
                    boundary_layers.append(layer_type(depth))

                return boundary_layers
        except Exception as e:
            logger.error(f"Error building boundary layers: {e}")
            raise

    def _build_cell(self, path: str) -> mp.Vector3:
        '''
        creates the simulations cell. Read directly from the config file
        
        Keywords: 
            path: location of the main config file
        
        Return:
            cell: mp.Vector3 of cell size 
        '''
        try:
            with open(path) as config_file:
                config = json.load(config_file)
                size_x = parse_expression(self.base_vars, config["cell"]["size_x"])
                size_y = parse_expression(self.base_vars, config["cell"]["size_y"])
                size_z = parse_expression(self.base_vars, config["cell"]["size_z"])

                return mp.Vector3(size_x, size_y, size_z)
        except Exception as e:
            logger.error(f"Error building simulation cell: {e}")
            raise

    def _build_flux(self, path : str) -> None:
        '''
        adds fluxes to the simulation domain if present
        '''
        try:
            with open(path) as fluxes:
                flux_config = json.load(fluxes)
                for flux in flux_config:
                    size = mp.Vector3(
                        parse_expression(self.base_vars, flux["size_x"]),
                        parse_expression(self.base_vars, flux["size_y"]),
                        parse_expression(self.base_vars, flux["size_z"]),
                    )
                    center = mp.Vector3(
                        parse_expression(self.base_vars, flux["center_x"]),
                        parse_expression(self.base_vars, flux["center_y"]),
                        parse_expression(self.base_vars, flux["center_z"]),
                    )

                    flux_region = mp.FluxRegion(center=center, size=size)
                    self.sim.add_flux(
                        self.base_vars["fcen"], self.base_vars["df"], self.base_vars["nfreq"], flux_region
                    )
        except Exception as e:
            logger.error(f"Error adding fluxes: {e}")

    def _load_ez_data(self, path : str) -> None:
        '''
        loads ez_data, if there is ez_data available in 
        
        Keywords:
            path: path to csv that stores the ez field after simulation across the simulation cell
            
        Return:
            ez field: 2d np.array of simulation domain with ez field strength at every postion
        '''
        return np.loadtxt(path, delimiter=",")

    def _append_to_config(self, params : dict) -> None:
        '''
        generic append function that allows to enlarge/edit the main config files by using a dictionary
        
        Keywords:
            params: dictionary of params that should be edited or changed
        '''
        try:
            # Load the existing configuration
            with open(self.config_file, 'r') as config_file:
                config = json.load(config_file)
                for key in params:
                    config[key] = params[key]
                    
            # Save the updated configuration back to the file
            with open(self.config_file, 'w') as config_file:
                json.dump(config, config_file, indent=4)
        except: 
            logger.error("File {} doesn't exist.".format(self.config_file))

    def _waveguide_ez_field(self, slice=None):
        '''
        gets the ez_data from the simulation and flat out returns it.
        '''
        self.ez_field = self.sim.get_array(center=mp.Vector3(), size=self.simulation_domain, component=mp.Ez).transpose()
        if slice == None:
            return self.ez_field
        else:
            return self.ez_field[slice]
    
    def _waveguide_ez_field_save(self, safe_file : str):
        '''
        after a simulation is done this saves the ez field to a csv with defined path
        
        Keywords:
            safe_file: location of the simulation ez field results
        '''
        ez_data = self._waveguide_ez_field()
        # write ez to a csv
        with open(safe_file,'w') as myfile:
            wr = csv.writer(myfile)
            wr.writerows(ez_data)
        # add ez file path to config directory
        self._append_to_config({"ez_file": safe_file})
        
#################################################################
# Multithreading support 
#################################################################

    def run(self, target, *args, **kwargs):
        """
        Start the simulation with a specified target function and its parameters.
    
        Keywords:
            target (callable): The function to execute in the thread.
            *args: Positional arguments for the target function.
            **kwargs: Keyword arguments for the target function.
        """
        if self._worker_thread and self._worker_thread.is_alive():
            raise RuntimeError("Simulation is already running.")
    
        self._stop_event.clear()
    
    
        def worker(target, *args, **kwargs):
            """
            Worker function that executes the target function with given parameters.
    
            Keywords:
                target (callable): The function to execute.
                *args: Positional arguments for the target function.
                **kwargs: Keyword arguments for the target function.
    
            Return:
                The result of the target function.
            """
            try:
                logger.info(f"Starting target function: {target.__name__}")
                result = target(*args, **kwargs)  # Execute the target function
                logger.info(f"Target function {target.__name__} completed successfully.")
                return result
            except Exception as e:
                logger.error(f"Error in target function {target.__name__}: {e}")
                raise
        # Start the worker thread
        
        self._worker_thread = threading.Thread(target=worker, args=(target, *args), kwargs=kwargs, daemon=False)
        self._worker_thread.start()

    def is_running(self) -> bool:
        """Check if the simulation is still running."""
        return self._worker_thread.is_alive() if self._worker_thread else False

    def stop(self):
        """Stop the simulation gracefully (if supported by Meep)."""
        if self._worker_thread and self._worker_thread.is_alive():
            self._stop_event.set()  # Signal the thread to stop
            self.sim.stop()  # Assuming Meep's stop method works cleanly

    def get_result(self):
        """Retrieve the simulation result if available."""
        if not self._result_queue.empty():
            return self._result_queue.get_nowait()
        else:
            raise RuntimeError("No result available yet. Simulation may still be running.")
      
      # Total simulation time

    def progress_callback(self, sim):
        self.progress_bar.update(1)

#################################################################
# External functions
#################################################################

    def view_sim_region(self) -> None:
        '''
        plots the simulation region to stdout. Useful to check config file.
        '''
        try:
            plt.figure(dpi=150)
            self.sim.plot2D()
            plt.show()
        except Exception as e:
            logger.error(f"Error visualizing simulation region: {e}")
 
    def get_waveguide_ez_field(self) -> np.array:
        '''
        external function to get the ez_field data from the simulation
        '''
        try:
            if self.ez_field.any() == None :
                # this will fail if ez_field has not been setup yet
                pass
        except:
            logger.warning("Simulation has not been run! Please run simulation")
            return [[0,0], [0,0]]
        else:
            return self.ez_field

        
# Derived class for ring resonator-specific functionality
class RingResonator(SimulationBase):
    '''
    Ring resonator simulation that offers the possibility to find resonances by harminv simulations and later offers
    full system simulation once a resonance for the specific simulation params have been found.
    '''
    def __init__(self, config_file : str) -> None:
        '''
        simple init that inherits from SimulationBase for simulation region setup
        then adds the ring resonator specific parameters of resonances to the system
        
        Keywords:
            config_file: location of the config file to set the simulation up
        '''
        super().__init__(config_file)
        self.resonance_frequency = -1
        self.resonance_wavelength = -1
  
    def find_resonances(self, resonance_conf_path : str):
        '''
        puts a harminv object within the ring resonator, resets the simulation (just in case the 
        simulation had some shenanigans done to it) and runs it
        
        Keywords:
            resonance_conf_path: output path to where to simulation details for the found resonance should be stored
        '''
        try:
            with open(self.geometry_path) as config_file:
                geometry_config = json.load(config_file)
                for geom in geometry_config:
                    if geom["type"] == "mp.Cylinder" and geom["radius"] == "radius": 
                        center_x = parse_expression(self.base_vars, geom["center_x"])
                        center_y = parse_expression(self.base_vars, geom["center_y"])
                        center_z = parse_expression(self.base_vars, geom["center_z"])
                        # Offset the position to a point near the ring's outer edge
                        y_offset = parse_expression(self.base_vars, geom["radius"]) + 0.5 * self.base_vars["ring_width"]
                        # Create a Harminv object for resonance analysis
                        h = mp.Harminv(
                            mp.Ez, 
                            mp.Vector3(center_x, center_y + y_offset, center_z),
                            fcen=self.base_vars["fcen"],
                            df=self.base_vars["df"]
                        )
                        self.sim.reset_meep()
                        # Run the simulation and collect resonance data
                        self.progress_bar = tqdm(total=self.runtime, desc="Simulation Progress")
                        self.sim.run(
                            mp.at_beginning(mp.output_epsilon),
                            mp.after_sources(h),
                            until_after_sources=self.runtime
                        )
                        self.progress_bar.close()
                        # Log the detected resonances
                        logger.info(f"Resonances detected: {h.modes}")
                        q = read_property_of_harminv(h, "Q")
                        resonances = read_property_of_harminv(h, "frequency")
                        self.resonance_frequency = self.get_resonance_from_q_factor(q, resonances)
                        self.resonance_wavelength = 1/self.resonance_frequency
                        self.write_simulation_config("conf_overcoupling.json", resonance_conf_path)
                        #self._append_to_config()
                        
        except Exception as e:
            logger.error(f"Error finding resonances: {e}")

    def get_resonance_from_q_factor(self, qFactor, resonances):
        '''
        find the resonance from the harminv q factor
        
        Keywords:
            qFactor: list of qFactors
            resonances: list of resonances
        Return:
            resonance: resonance frequency of the ring if it exists
        '''
        # Normalize Q factors
        max_q = max(qFactor)
        normalized_q = [q / max_q for q in qFactor]


        # plt.scatter(resonances, normalized_q, label="Q factors", color="blue")
        # plt.xlabel("Frequency")
        # plt.ylabel("Q Factor")
        # plt.show()
        # Filter modes with significant Q factors (e.g., >50% of max Q)
        threshold = 0.5
        significant_indices = [i for i, q in enumerate(normalized_q) if q > threshold]
        significant_resonances = [resonances[i] for i in significant_indices]
        significant_q = [qFactor[i] for i in significant_indices]

        # Fit Gaussian to the significant Q factors
        if len(significant_resonances) > 3:  # Need at least 3 points for a meaningful fit
            try:
                # Perform curve fitting
                popt, _ = curve_fit(gaussian, significant_resonances, significant_q, p0=[1, resonances[np.argmax(qFactor)], 0.1])
                fitted_amp, fitted_mean, fitted_std = popt

                # Plot the Q factors and Gaussian fit
                x_fit = np.linspace(min(significant_resonances), max(significant_resonances), 500)
                y_fit = gaussian(x_fit, *popt)

                #plt.scatter(significant_resonances, significant_q, label="Significant Q factors", color="blue")
                #plt.plot(x_fit, y_fit, label="Gaussian Fit", color="red")
                #plt.xlabel("Frequency")
                #plt.ylabel("Q Factor")
                #plt.legend()
                #plt.show()

                # Choose the frequency corresponding to the Gaussian mean
                chosen_frequency = fitted_mean
                print(f"Chosen frequency based on Gaussian fit: {chosen_frequency}")
                return chosen_frequency
            except RuntimeError as e:
                print(f"Gaussian fit failed: {e}")
                chosen_frequency = resonances[np.argmax(qFactor)]
                print(f"Defaulting to frequency with highest Q factor: {chosen_frequency}")
                return chosen_frequency
        else:
            # Default to highest Q factor if Gaussian fit is not viable
            chosen_frequency = resonances[np.argmax(qFactor)]
            print(f"Chosen frequency with highest Q factor: {chosen_frequency}")
            return chosen_frequency

    def write_simulation_config(self, base_config_file : str, new_file_name: str):
        """
        Updates the resonance wavelength in the 'base_vars' of the configuration file.

        Keywords:
            config_file: Path to the configuration file.
            new_wavelength: The new resonance wavelength to be updated.
        """
        try:
            # Load the existing configuration
            with open(base_config_file, 'r') as f:
                config = json.load(f)

            # Update the resonance wavelength in base_vars
            if "base_vars" in config:
                config["base_vars"]["wavelength"] = self.resonance_wavelength
                # this should not be constant but idk what a dynamic value might be, so it stays
                config["base_vars"]["df"] = 0.02
            else:
                raise KeyError("The 'base_vars' section is missing in the configuration file.")
            if "source_file" in config:
                config["source_file"] = "simulation_sources.json"
            else:
                raise KeyError("The 'source_file' section is missing in the configuration file.")
            # Recompute derived_base_vars if it exists
            if "derived_base_vars" in config and "fcen" in config["derived_base_vars"]:
                config["derived_base_vars"]["fcen"] = f"1/{self.resonance_wavelength}"

            # Save the updated configuration back to the file
            with open(new_file_name, 'w') as f:
                json.dump(config, f, indent=4)

            print(f"Configuration updated successfully. New wavelength: {self.resonance_wavelength}")

        except FileNotFoundError:
            print(f"Error: File {base_config_file} not found.")
        except KeyError as e:
            print(f"Error: Missing key in configuration: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

    def run_resonance(self, save_ez=False):
        """
        once resonances are found, this function runs the simulation until the self.runtime is reached
        
        Keywords:
            save_ez: boolean of whether you want to save ez or just discard it
        """
        # TODO: if the resonance has been run already: ask if it should be rerun
        
        # check if the source is the correct 
        try:
            with open(self.sources_path) as f:
                sources = json.load(f)
                # allows to run the simulation ONLY if we're loading the correct sources - else aborts
                valid_simulation_sources = [True for source in sources if source["usage"] == "simulation"]
                if (len(valid_simulation_sources) == len(sources) and (len(sources) != 0)):
                    self.progress_bar = tqdm(total=self.runtime, desc="Simulation Progress")
                    self.sim.run(mp.at_beginning(mp.output_epsilon),
                                 mp.to_appended("{}".format(self.config_file.split(".")[0]), mp.at_every(1, lambda sim: self.progress_callback(sim))),
                                 until=self.runtime)
                    self.progress_bar.close()
                    # save ez to a file
                    safe_file = ""
                    if save_ez == True:
                        for elem in self.config_file.split("."):
                            if elem != self.config_file.split(".")[-1]:
                                safe_file += elem
                                #safe_file += "."
                        safe_file = "ez_" + safe_file + "csv"
                        self._waveguide_ez_field_save(safe_file)
                else:
                    logger.error("Check {}! Not all source are having the right usage".format(self.sources_path))
        except:
            raise FileExistsError("File {} does not exist".format(self.sources_path))

    def plot_field_result(self, slice: tuple) -> None:
        '''
        plot the field results 
        
        Keywords:
            slice: tuple of a y axis section that should be in focus while plotting
        '''
        plt.imshow(self.ez_field[slice[0]:slice[1]])
        plt.show()

    def plot_coupling(self, plot_params_guess=[0.08,0.3, 5], slice=90):
        '''
        rudimentary fitting function that should take the wavefront along the linear waveguide and then fits a sinusoidal
        function to the interferred part of the function
        
        Keywords:
            plot_params_guess: guess for the fitting parameters (default works well for the set example )
            slice: 1d array element of the row in which the linear waveguide sits within the simulation domain
        '''
        ez_data = self.sim.get_array(center=mp.Vector3(), size=self.simulation_domain, component=mp.Ez)
        plt.plot(ez_data.transpose()[slice],  "--", alpha=0.5, label="simulation")
        xdata_start= 250
        xdata_end = 370
        ydata = ez_data.transpose()[90][xdata_start:xdata_end]
        xdata = np.arange(xdata_start, xdata_end)
        guess = plot_params_guess
        popt, pcov = curve_fit(sinus, xdata, ydata, p0=guess)
        fit_x = np.arange(0, len(ez_data.transpose()[90]))
        plt.plot(fit_x, sinus(fit_x, *popt), label="extrapolation")
        plt.plot(xdata, ydata, label='fit range')
        plt.legend()
        plt.xlabel("distance along linear waveguide [arb. u.]")
        plt.ylabel("intensity of e field [arb. u.]")
        plt.plot()
        plt.grid()

def sinus(x, a, b, c):
    return a * np.sin(b * x + c)

# Define a Gaussian function for fitting
def gaussian(x, amp, mean, std_dev):
    return amp * np.exp(-((x - mean) ** 2) / (2 * std_dev ** 2))

def read_property_of_harminv(harminv : mp.Harminv, key): 
    property = []
    for mode in harminv.modes:
        # frequency
        if key == "frequency":
            property.append(mode.freq) 
        # imag. freq.
        elif key == "imag. freq.":
            property.append(mode.amp.imag) 
        # Q
        elif key == "Q":
            property.append(mode.Q) 
        # |amp|
        elif key == "abs amplitude":
            property.append(np.sqrt(mode.amp.real**2))
        # amplitude
        elif key == "amplitude":
            property.append(np.sqrt(mode.amp.real))
        # error
        elif key == "error":
            property.append(mode.err)
        else: 
            raise Exception("Unkown Key")

    return property
