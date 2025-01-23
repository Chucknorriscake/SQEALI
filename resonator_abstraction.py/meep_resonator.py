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

# Configure logging for debugging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Base class for handling configuration and setup
class SimulationBase:
    

    def __init__(self, config_file=None):
        self.base_vars = {}
        self.sources = []
        self.geometry = []
        self.boundary_layers = []
        self.simulation_domain = mp.Vector3(10, 10, 0)
        self.resolution = 16
        self.runtime = 100
        self.sim = None
        self._worker_thread = None
        mp.divide_parallel_processes(1)
        # shut up
        mp.verbosity.meep = 0
        self._result_queue = queue.Queue()
        self._stop_event = threading.Event()

        
        if config_file:
            self.sim, self.sources_path ,self.flux_path,self.geometry_path = self._load_conf(config_file)

    def _load_conf(self, path: str) -> mp.Simulation:
        try:
            
            with open(path) as config_file:
                config = json.load(config_file)
                sources_path = config["source_file"]
                flux_path = config["flux_file"]
                geometry_path = config["geometry_file"]
            
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

    def _set_simulation_params(self):
        self.resolution = self.base_vars.get("resolution", 16)
        self.runtime = self.base_vars.get("runtime", 100)

    def _build_sources(self, path: str) -> list:
        try:
            with open(path) as f:
                sources_config = json.load(f)
                print(sources_config)
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

    def view_sim_region(self) -> None:
        try:
            plt.figure(dpi=150)
            self.sim.plot2D()
            plt.show()
        except Exception as e:
            logger.error(f"Error visualizing simulation region: {e}")

    def run(self, target, *args, **kwargs):
        """
        Start the simulation with a specified target function and its parameters.
    
        Args:
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
    
            Args:
                target (callable): The function to execute.
                *args: Positional arguments for the target function.
                **kwargs: Keyword arguments for the target function.
    
            Returns:
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
      

# Derived class for ring resonator-specific functionality
class RingResonator(SimulationBase):
    def __init__(self, config_file=None):
        super().__init__(config_file)
        self.resonance_frequency = -1
        self.resonance_wavelength = -1
        self.add_fluxes()

    def add_fluxes(self):
        try:
            with open(self.flux_path) as fluxes:
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

    def find_resonances(self, resonance_conf_path : str):
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
                        self.sim.run(
                            mp.at_beginning(mp.output_epsilon),
                            mp.after_sources(h),
                            until_after_sources=self.runtime
                        )
                        
                        # Log the detected resonances
                        logger.info(f"Resonances detected: {h.modes}")
                        q = read_property_of_harminv(h, "Q")
                        resonances = read_property_of_harminv(h, "frequency")
                        self.resonance_frequency = self.get_resonance_from_q_factor(q, resonances)
                        self.resonance_wavelength = 1/self.resonance_frequency
                        self.write_simulation_config("conf_overcoupling.json", resonance_conf_path)
                        
        except Exception as e:
            logger.error(f"Error finding resonances: {e}")

    def get_resonance_from_q_factor(self, qFactor, resonances):

        # Normalize Q factors
        max_q = max(qFactor)
        normalized_q = [q / max_q for q in qFactor]


        plt.scatter(resonances, normalized_q, label="Q factors", color="blue")
        plt.xlabel("Frequency")
        plt.ylabel("Q Factor")
        plt.show()
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

                plt.scatter(significant_resonances, significant_q, label="Significant Q factors", color="blue")
                plt.plot(x_fit, y_fit, label="Gaussian Fit", color="red")
                plt.xlabel("Frequency")
                plt.ylabel("Q Factor")
                plt.legend()
                plt.show()

                # Choose the frequency corresponding to the Gaussian mean
                chosen_frequency = fitted_mean
                return chosen_frequency
                print(f"Chosen frequency based on Gaussian fit: {chosen_frequency}")
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

        Args:
            config_file (str): Path to the configuration file.
            new_wavelength (float): The new resonance wavelength to be updated.
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

    def run_resonance(self):
        # check if the source is the correct 
        try:
            with open(self.sources_path) as f:
                sources = json.load(f)
                print(sources)
                valid_simulation_sources = [True for source in sources if source["useage"] == "simulation"]
                if (len(valid_simulation_sources) == len(sources) and (len(sources) != 0)):
                    self.sim.run(mp.at_beginning(mp.output_epsilon),
                                 mp.to_appended("{}".format(self.config_file.split(".")[0]), mp.at_every(1, mp.output_efield_z)),
                                 until=self.runtime)
                    pass
                else:
                    raise ValueError("Check {}! Not all source are having the right usage".format(self.sources_path))
        except:
            raise FileExistsError("File {} does not exist".format(self.sources_path))

    def plot_field_result(self):
        ez_data = self.sim.get_array(center=mp.Vector3(), size=self.simulation_domain, component=mp.Ez)
        plt.imshow(ez_data.transpose()[70:200])
        plt.show()

    def plot_coupling(self, plot_params_guess):
        ez_data = self.sim.get_array(center=mp.Vector3(), size=self.simulation_domain, component=mp.Ez)
        plt.plot(ez_data.transpose()[90],  "--", alpha=0.5, label="simulation")


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

    def waveguide_ez_field(self, slice=None):
        ez_data = self.sim.get_array(center=mp.Vector3(), size=self.simulation_domain, component=mp.Ez)
        if slice == None:
            return ez_data.transpose()
        else:
            return ez_data.transpose()[slice]
        
        
        
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
