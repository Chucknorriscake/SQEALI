import numpy as np
import matplotlib.pyplot as plt
import meep as mp
import json
import logging
import threading
import os
import queue
from IPython.display import Video
from scipy.optimize import curve_fit
from scipy.signal import hilbert, argrelextrema
from helpers import parse_expression
from tqdm import tqdm
import csv 
# Configure logging for debugging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Base class for handling configuration and setup
class SimulationBase:
    """
    Meep wrapper that allows building simulations by loading a JSON config file and 
    setting up the simulation domain accordingly.
    """
    def __init__(self, config_file : str) -> None:
        """
        Initializes the simulation environment using a given JSON config file.

        Args:
            config_file (str): Path to the configuration file containing simulation parameters.
        """
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
        self.ez_loaded = False
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
        """
        Loads the configuration from a JSON file and sets up the simulation.

        Args:
            path (str): Path to the configuration file.

        Returns:
            mp.Simulation: Configured Meep simulation object.
        """
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
                    self.ez_loaded = True
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
        """
        Loads basic and derived variables from the configuration.

        Args:
            path (str): Path to the configuration file.

        Returns:
            dict: Dictionary containing basic simulation variables.
        """
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
        """
        Sets resolution and runtime parameters for the simulation.
        """
        self.resolution = self.base_vars.get("resolution", 16)
        self.runtime = self.base_vars.get("runtime", 100)

    def _build_sources(self, path: str) -> list:
        """
        Builds sources from a JSON file.

        Args:
            path (str): Path to the source configuration file.

        Returns:
            list: List of Meep source objects.
        """
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
        """
        Builds geometry objects from a configuration file.

        Args:
            path (str): Path to the geometry configuration file.

        Returns:
            list: List of Meep geometry objects (Blocks or Cylinders).
        """
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
        """
        Builds boundary layers from the main configuration file.

        Args:
            path (str): Path to the configuration file.

        Returns:
            list: List of Meep boundary layer objects.
        """
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
        """
        Constructs the simulation cell dimensions from config.

        Args:
            path (str): Path to the configuration file.

        Returns:
            mp.Vector3: Vector describing the simulation cell size.
        """
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
        """
        Adds flux monitors to the simulation from a configuration file.

        Args:
            path (str): Path to the flux configuration file.
        """
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
        """
        Loads electric field (Ez) data from a CSV file.

        Args:
            path (str): Path to the CSV file.

        Returns:
            np.ndarray: 2D array of Ez field values.
        """
        return np.loadtxt(path, delimiter=",")

    def _append_to_config(self, params : dict) -> None:
        """
        Updates the configuration file with new parameters.

        Args:
            params (dict): Dictionary of parameters to append/update in the config.
        """
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
        """
        Retrieves the Ez field from the simulation domain.

        Args:
            slice (optional): Slice of the Ez field to return.

        Returns:
            np.ndarray: 2D array or sliced view of the Ez field.
        """
        self.ez_field = self.sim.get_array(center=mp.Vector3(), size=self.simulation_domain, component=mp.Ez).transpose()
        if slice == None:
            return self.ez_field
        else:
            return self.ez_field[slice]
    
    def _waveguide_ez_field_save(self, safe_file : str):
        """
        Saves the Ez field data to a CSV file.

        Args:
            safe_file (str): Path to the file where Ez data will be saved.
        """
        ez_data = self._waveguide_ez_field()
        # write ez to a csv
        with open(safe_file,'w') as myfile:
            wr = csv.writer(myfile)
            wr.writerows(ez_data)
        # add ez file path to config directory
        self._append_to_config({"ez_file": safe_file})
        
    def run_simulation(self, save_ez=False):
        """
        Runs the simulation using the configured parameters.

        Args:
            save_ez (bool, optional): Whether to save the Ez field data. Defaults to False.
        """
        try:
            with open(self.sources_path) as f:
                sources = json.load(f)
                # allows to run the simulation ONLY if we're loading the correct sources - else aborts
                valid_simulation_sources = [True for source in sources if source["usage"] == "simulation"]
                if (len(valid_simulation_sources) == len(sources) and (len(sources) != 0)):
                    clear_cmd()
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
                                safe_file += "."
                        safe_file = "ez_" + safe_file + "csv"
                        self._waveguide_ez_field_save(safe_file)
                else:
                    logger.error("Check {}! Not all source are having the right usage".format(self.sources_path))
        except:
            raise FileExistsError("File {} does not exist".format(self.sources_path))

#################################################################
# Multithreading support 
#################################################################

    def run(self, target, *args, **kwargs):
        """
        Runs a target function in a separate thread.

        Args:
            target (callable): Function to execute.
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
        """
        Checks if the worker thread is still running.

        Returns:
            bool: True if the thread is active, False otherwise.
        """
        return self._worker_thread.is_alive() if self._worker_thread else False

    def stop(self):
        """Stop the simulation gracefully (if supported by Meep)."""
        if self._worker_thread and self._worker_thread.is_alive():
            self._stop_event.set()  # Signal the thread to stop
            self.sim.stop()  # Assuming Meep's stop method works cleanly

    def get_result(self):
        """
        Retrieves the result from the result queue.

        Returns:
            any: Result from the simulation.

        Raises:
            RuntimeError: If no result is available yet.
        """
        if not self._result_queue.empty():
            return self._result_queue.get_nowait()
        else:
            raise RuntimeError("No result available yet. Simulation may still be running.")
      
      # Total simulation time

    def progress_callback(self, sim):
        """
        Updates the progress bar during simulation.

        Args:
            sim (mp.Simulation): The Meep simulation object.
        """
        self.progress_bar.update(1)

#################################################################
# External functions
#################################################################

    def view_sim_region(self, dpi=150) -> None:
        """
        Plots the simulation region.

        Args:
            dpi (int, optional): DPI resolution for the plot. Defaults to 150.
        """
        try:
            plt.figure(dpi=dpi)
            self.sim.plot2D()
            plt.show()
        except Exception as e:
            logger.error(f"Error visualizing simulation region: {e}")
 
    def get_waveguide_ez_field(self) -> np.array:
        """
        Gets the Ez field after simulation.

        Returns:
            np.ndarray: 2D array of the Ez field.
        """
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
    """
    Ring resonator simulation that offers the possibility to find resonances via Harminv
    simulations and later supports full system simulation after a resonance is identified.
    """
    def __init__(self, config_file : str) -> None:
        """
        Initializes the ring resonator simulation by inheriting from SimulationBase.

        Args:
            config_file (str): Path to the configuration file.
        """
        super().__init__(config_file)
        self.resonance_frequency = -1
        self.resonance_wavelength = -1
  
    def find_resonances(self, resonance_conf_path : str):
        """
        Runs a Harminv-based simulation to find resonances in the ring resonator. 
        The resonance frequency and wavelength are extracted and saved.

        Args:
            resonance_conf_path (str): Path to save the updated simulation configuration file.
        """
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
                        clear_cmd()
                        self.sim.run(
                            mp.at_beginning(mp.output_epsilon),
                            mp.after_sources(h),
                            mp.at_every(1, lambda sim: self.progress_callback(sim)),
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
        """
        Determines the resonance frequency based on the Q factors returned by Harminv.

        Args:
            qFactor (list): List of Q factor values.
            resonances (list): Corresponding list of resonance frequencies.

        Returns:
            float: The selected resonance frequency.
        """
        # Normalize Q factors
        max_q = max(qFactor)
        normalized_q = [q / max_q for q in qFactor]

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
        Updates the resonance wavelength in the 'base_vars' of the given configuration file
        and writes the modified configuration to a new file.

        Args:
            base_config_file (str): Path to the base configuration file.
            new_file_name (str): Path to save the updated configuration.
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
        Runs a full simulation using the previously determined resonance frequency.

        Args:
            save_ez (bool, optional): Whether to save the Ez field data. Defaults to False.
        """
        # TODO: if the resonance has been run already: ask if it should be rerun
        
        # check if the source is the correct 
        try:
            with open(self.sources_path) as f:
                sources = json.load(f)
                # allows to run the simulation ONLY if we're loading the correct sources - else aborts
                valid_simulation_sources = [True for source in sources if source["usage"] == "simulation"]
                if (len(valid_simulation_sources) == len(sources) and (len(sources) != 0)):
                    clear_cmd()
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
                                safe_file += "."
                        safe_file = "ez_" + safe_file + "csv"
                        self._waveguide_ez_field_save(safe_file)
                else:
                    logger.error("Check {}! Not all source are having the right usage".format(self.sources_path))
        except:
            raise FileExistsError("File {} does not exist".format(self.sources_path))

    def plot_field_result(self, slice: tuple) -> None:
        """
        Plots the Ez field result using a vertical slice of the simulation domain.

        Args:
            slice (tuple): A tuple specifying the start and end indices along the y-axis.
        """
        plt.imshow(self.ez_field[slice[0]:slice[1]])
        plt.show()

    def plot_coupling(self, plot_params_guess=[0.08,0.3, 5], slice=90):
        """
        Fits a sinusoidal curve to the wavefront along the linear waveguide and plots the results.

        Args:
            plot_params_guess (list, optional): Initial guess for the sinusoidal fitting parameters.
            slice (int, optional): Y-axis index corresponding to the waveguide cross-section.
        """
        ez_data = self.get_waveguide_ez_field() #self.sim.get_array(center=mp.Vector3(), size=self.simulation_domain, component=mp.Ez)
        plt.plot(ez_data.transpose()[slice],  "--", alpha=0.5, label="simulation")
        xdata_start= 250
        xdata_end = 370
        ydata = ez_data.transpose()[90]#[xdata_start:xdata_end]
        xdata = np.arange(0, len(ydata))#np.arange(xdata_start, xdata_end)
        guess = plot_params_guess
        try:
            popt, pcov = curve_fit(sinus, xdata, ydata, p0=guess)
            fit_x = np.arange(0, len(ez_data.transpose()[90]))
            #plt.plot(fit_x, sinus(fit_x, *popt), label="extrapolation")
            #plt.plot(xdata, ydata, label='fit range')
        except:
            logger.info("No fit found")
            #plt.plot(xdata, ydata, label='fit range')
        plt.legend()
        plt.xlabel("distance along linear waveguide [arb. u.]")
        plt.ylabel("intensity of e field [arb. u.]")
        plt.plot()
        plt.grid()
        


class LinearWaveguides(SimulationBase):
    """
    A class for simulating and analyzing coupling behavior between parallel linear waveguides.
    """
    def __init__(self, config_file):
        """
        Initializes the LinearWaveguides simulation.

        Args:
            config_file (str): Path to the configuration file used for setting up the simulation.
        """
        super().__init__(config_file)
        
    
    def calc_coupling_length(self) -> float:
        """
        Calculates the coupling length between waveguides based on the minima of the field envelope.

        This method uses the Hilbert transform to extract the amplitude envelope of the Ez field
        and identifies local minima that fall below a certain threshold, defined as a percentage of the
        maximum envelope value. The average distance between these filtered minima is then used to
        compute the coupling length.

        Returns:
            float: The calculated coupling length in simulation units.
        """
        filter = []
        threshold_percent = 0.1
        waveguide_center = int(np.round((self.base_vars["sy"] + self.base_vars["waveguide_distance"])/2*self.base_vars["resolution"],0))
        signal = self.ez_field[waveguide_center]
        # Compute analytic signal and amplitude envelope
        analytic_signal = hilbert(signal)
        amplitude_envelope = np.abs(analytic_signal)

        # Find all local minima
        min_indices = argrelextrema(amplitude_envelope, np.less)[0]


        with open(self.config_file, 'r') as config_file:
            config = json.load(config_file)
            threshold_percent = config["coupling_fit_threshold"]

        with open(self.config_file, 'r') as config_file:
            config = json.load(config_file)
            filter.append(config["coupling_fit_first"])
            filter.append(config["coupling_fit_last"])
        
        # Threshold: keep only minima below 10% of the max envelope value
        threshold = threshold_percent * np.max(amplitude_envelope)
        filtered_min_indices = min_indices[amplitude_envelope[min_indices] < threshold][filter[0]:filter[1]]
        coupling_length = np.average(np.diff(filtered_min_indices))/self.base_vars["resolution"]
        return coupling_length
        
    def plot_coupling_length(self, threshold_percent : float = None ,filter: list = []):
        """
        Plots the Ez field and its amplitude envelope, highlighting the detected local minima used 
        in coupling length estimation.

        This is useful for visual debugging and validation of coupling length calculations.
        
        Changes in threshold_percent and filter will automatically be updated to the config file

        Args:
            threshold_percent (float, optional): Custom threshold as a percentage of max envelope for minima selection.
                If None, the value is read from the configuration file.
            filter (list, optional): List containing two indices `[start, end]` to slice the set of minima. 
                If not provided, values are taken from the configuration.
        """
        waveguide_center = int(np.round((self.base_vars["sy"] + self.base_vars["waveguide_distance"])/2*self.base_vars["resolution"],0))
        signal = self.ez_field[waveguide_center]
        # Compute analytic signal and amplitude envelope
        analytic_signal = hilbert(signal)
        amplitude_envelope = np.abs(analytic_signal)

        # Find all local minima
        min_indices = argrelextrema(amplitude_envelope, np.less)[0]


        # filtering
        if threshold_percent == None:
            with open(self.config_file, 'r') as config_file:
                config = json.load(config_file)
                threshold_percent = config["coupling_fit_threshold"]
                print(threshold_percent)
        else: 
            new_params = {"coupling_fit_threshold" : threshold_percent}
            self._append_to_config(new_params)
            
        if len(filter) == 0:
            with open(self.config_file, 'r') as config_file:
                config = json.load(config_file)
                filter.append(config["coupling_fit_first"])
                filter.append(config["coupling_fit_last"])
        else:
            new_params = {
                "coupling_fit_first":filter[0],
                "coupling_fit_last":filter[1],
            }
            self._append_to_config(new_params)

        # Threshold: keep only minima below 10% of the max envelope value
        threshold = threshold_percent * np.max(amplitude_envelope)
        filtered_min_indices = min_indices[amplitude_envelope[min_indices] < threshold][filter[0]:filter[1]]

        # Plot
        plt.plot(amplitude_envelope, label="Envelope")
        plt.plot(signal, label="Original Signal")
        plt.scatter(filtered_min_indices, amplitude_envelope[filtered_min_indices], color='red', label='Minima < 10% max')
        plt.legend()
        plt.title("Envelope Minima Below 10% Threshold")
        print(np.diff(filtered_min_indices))
        plt.show()
        
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


def clear_cmd():
    os.system('clear')