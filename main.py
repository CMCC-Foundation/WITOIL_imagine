# ------------------------------------------------
# MEDSLIK-II oil spill fate and transport model
# ------------------------------------------------
import sys, os
import shutil
import logging
import warnings
import subprocess
import logging
import importlib.util
import numpy as np
import pandas as pd
from datetime import datetime
from argparse import ArgumentParser
from glob import glob

# Import medslik modules
from src.utils.config import Config
from src.utils.utils import Utils
from src.download.download_era5_parser import *
from src.download.download_copernicus_parser import *
from src.preprocessing import PreProcessing
from src.postprocessing import PostProcessing
from src.plot import MedslikIIPlot
from src.utils.read_oil_data import read_oilbase
from src.utils.utils import Utils

# Import pyngl if present
package_spec = importlib.util.find_spec("Ngl")
if package_spec is not None:
    import Ngl

    _has_pyngl = True
else:
    _has_pyngl = False

# Create a logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Create a file handler with overwrite mode ('w')
file_handler = logging.FileHandler("medslik_run.log", mode="w")
file_handler.setLevel(logging.DEBUG)

# Create a formatter and set it for the handler
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
file_handler.setFormatter(formatter)

# Add the file handler to the logger
logger.addHandler(file_handler)


class MedslikII:
    """
    This class embeds the MAIN code of medslik-II software.
    """

    def __init__(self, config: dict) -> None:
        """
        Class constructor given config file path.
        """
        self.config = config
        self.medslik_directory = os.path.dirname(os.path.abspath(__file__))
        # Create experiment directories
        self.root_directory = os.path.join(
            self.config["simulation"]["experiment_path"], config["simulation"]["name"]
        )
        os.makedirs(self.root_directory, exist_ok=True)
        self.out_directory = os.path.join(self.root_directory, "out_files")
        os.makedirs(self.out_directory, exist_ok=True)
        self.out_figures = os.path.join(self.out_directory, "figures")
        os.makedirs(self.out_figures, exist_ok=True)
        self.xp_directory = os.path.join(self.root_directory, "xp_files")
        os.makedirs(self.xp_directory, exist_ok=True)
        spill_lat = np.array(self.config["simulation"]["spill_lat"])
        self.n_spill_points = np.shape(spill_lat)[0]
        if config["input_files"]["set_domain"] == True:
            lat_min, lat_max = (
                config["input_files"]["lat"][0],
                config["input_files"]["lat"][1],
            )
            lon_min, lon_max = (
                config["input_files"]["lon"][0],
                config["input_files"]["lon"][1],
            )
        else:
            latitude = config["simulation"]["spill_lat"][0]
            longitude = config["simulation"]["spill_lon"][0]

            lat_min, lat_max = (
                latitude - config["input_files"]["delta"][0],
                latitude + config["input_files"]["delta"][0],
            )
            lon_min, lon_max = (
                longitude - config["input_files"]["delta"][0],
                longitude + config["input_files"]["delta"][0],
            )
        self.lat_min, self.lat_max = lat_min, lat_max
        self.lon_min, self.lon_max = lon_min, lon_max
        self.apply_aging_effects()
        self.initial_checking()

    def apply_aging_effects(self) -> None:
        """
        Consider aging when setting
        the simulation length and the start datetime.
        """
        oilspill = self.config["simulation"]
        age = oilspill["slick_age"]
        if self.n_spill_points > 1:
            oilspill["sim_length"] += age
            oilspill["start_datetime"] -= pd.Timedelta(hours=age)
            oilspill["spill_duration"] = 0
        self.config["simulation"] = oilspill

    def initial_checking(self):
        """
        Check if any issue might derive from configuration.
        """
        # checking if the coordinates are on land
        lat = self.config["simulation"]["spill_lat"]
        lon = self.config["simulation"]["spill_lon"]
        sea = Utils.check_land(lon, lat)
        if sea == 0:
            raise ValueError(
                "Your coordinates lie within land. Please check your values again"
            )
        # checking dates
        dt = Utils.validate_date(self.config["simulation"]["start_datetime"])
        logger.info("No major issues found on dates and oil spill coordinates")

        # checking it starting from area spill
        if self.config["input_files"]["shapefile"]["shape_path"]:
            shapefile_path = self.config["input_files"]["shapefile"]["shape_path"]
            if os.path.exists(shapefile_path):
                logger.info(
                    f"Simulation initial conditions area spill are provided on \
                        {self.config['input_files']['shapefile']['shape_path']}. \
                        Spill rate from config files will not be considered"
                )
                volume = Utils.oil_volume_shapefile(self.config)

                # Correcting volume on the config object
                self.config["simulation"]["spill_rate"] = volume

    @staticmethod
    def data_download_medslik(
        config: dict, domain: list[float], root_directory: str
    ) -> None:
        """
        Download METOCE datasets.
        """
        lon_min, lon_max, lat_min, lat_max = domain
        copernicus_user = config["download"]["copernicus_user"]
        copernicus_pass = config["download"]["copernicus_password"]
        date = pd.to_datetime(config["simulation"]["start_datetime"])
        identifier = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
        inidate = date - pd.Timedelta(hours=1)
        enddate = date + pd.Timedelta(hours=config["simulation"]["sim_length"] + 24)
        if (
            30.37 < float(lat_min) < float(lat_max) < 45.7
            and -17.25 < float(lon_min) < float(lon_max) < 35.9
        ):
            down = "local"
        else:
            down = "global"
        if config["download"]["download_curr"]:
            output_path = "data/COPERNICUS/"
            output_name = output_path + "Copernicus{}_{}_{}_mdk.nc".format(
                "{}", identifier, config["simulation"]["name"]
            )
            logger.info("Downloading CMEMS currents")
            download_copernicus(
                lat_min,
                lat_max,
                lon_min,
                lon_max,
                0,
                120,
                inidate,
                enddate,
                down,
                output_path=output_path,
                output_name=output_name,
                user=copernicus_user,
                password=copernicus_pass,
            )
            subprocess.run(
                [
                    f'cp {output_path}*{identifier}*{config["simulation"]["name"]}*.nc {root_directory}/oce_files/'
                ],
                shell=True,
            )
            subprocess.run(
                [f'rm {output_path}*{identifier}*{config["simulation"]["name"]}*.nc'],
                shell=True,
            )
        if config["download"]["download_wind"]:
            output_path = "data/ERA5/"
            output_name = output_path + "era5_winds10_{}_{}_mdk.nc".format(
                identifier, config["simulation"]["name"]
            )
            logger.info("Downloading ERA5 reanalysis winds")
            get_era5(
                lon_min,
                lon_max,
                lat_min,
                lat_max,
                inidate,
                enddate,
                output_path=output_path,
                output_name=output_name,
            )
            process_era5(output_path=output_path, output_name=output_name)
            subprocess.run(
                [
                    f'cp {output_path}*{identifier}*{config["simulation"]["name"]}*.nc {root_directory}/met_files/'
                ],
                shell=True,
            )
            subprocess.run(
                [f'rm {output_path}*{identifier}*{config["simulation"]["name"]}*.nc'],
                shell=True,
            )

    @staticmethod
    def run_preproc(config: dict, exp_folder: str, lon_min, lon_max, lat_min, lat_max):
        """
        Run preprocessing.
        """
        domain = [lon_min, lon_max, lat_min, lat_max]
        preproc = PreProcessing(config=config, exp_folder=exp_folder, domain=domain)
        # Create folders
        preproc.create_directories()
        # download data if needed
        if config["download"]["download_data"] == True:
            MedslikII.data_download_medslik(config, domain, exp_folder)
        if config["run_options"]["preprocessing"]:
            if config["run_options"]["preprocessing_metoce"]:
                oce_path = config["input_files"]["metoce"]["oce_data_path"]
                met_path = config["input_files"]["metoce"]["met_data_path"]
                if oce_path == "":
                    oce_path = None
                if met_path == "":
                    met_path = None
                preproc.process_currents(oce_path=oce_path)
                preproc.process_winds(met_path=met_path)
            if config["run_options"]["preprocessing_dtm"]:
                preproc.common_grid()
                preproc.process_bathymetry(
                    config["input_files"]["dtm"]["bathymetry_path"]
                )
                preproc.process_coastline(
                    config["input_files"]["dtm"]["coastline_path"]
                )
            if config["simulation"]["area_spill"]:
                if config["input_files"]["shapefile"]["shape_path"]:
                    shapefile_path = config["input_files"]["shapefile"]["shape_path"]
                    if os.path.exists(shapefile_path):
                        # using an area spill identified on a shapefile to generate initial conditions
                        preproc.process_initial_shapefile()
                else:
                    # passing the area vertex to generate the polygon
                    preproc.write_initial_input(config["simulation"]["area_vertex"])

    @staticmethod
    def toml_to_parameters(path_to_toml: str, txt_file: str) -> None:
        """
        Write parameters.txt from parameters.toml file.
        """
        toml_file = Config(path_to_toml).config_dict
        f = open(txt_file, mode="w")
        for key, keyval in toml_file.items():
            f.write(key + "\n")
            for key, val in keyval.items():
                if type(val) == list:
                    string_val = " ".join(map(str, val))
                else:
                    string_val = str(val)
                if string_val == "True":
                    string_val = "1"
                if string_val == "False":
                    string_val = "0"
                f.write(string_val + "\n")
        f.close()

    @staticmethod
    def config_to_input(
        config: dict, oil_file_path: str, path_to_inp_file: str
    ) -> None:
        """
        Convert config dictionary into medslik.inp file.
        """
        current_path = os.path.dirname(os.path.abspath(__file__))
        template_file = os.path.join(current_path, "src", "model", "oilspill.inp")
        simulation = config["simulation"]

        # Checks if spill starts from area spill or point source
        if simulation["area_spill"]:
            isat = 1
        else:
            isat = 0

        replace_dict = {
            "restart_hr": "0",
            "restart_min": "0",
            "day": simulation["start_datetime"].day,
            "month": simulation["start_datetime"].month,
            "year": simulation["start_datetime"].year,
            "hour": simulation["start_datetime"].hour,
            "minutes": f"{simulation['start_datetime'].minute:02d}",
            "iage": simulation["slick_age"],
            "isat": isat,
            "spl_dur": int(simulation["spill_duration"][0]),
            "sim_length": int(simulation["sim_length"]),
            "spl_lat": simulation["spill_lat"],
            "spl_lon": simulation["spill_lon"],
            "splrate": simulation["spill_rate"],
        }
        with open(template_file, "r") as file:
            file_contents = file.read()
        for key, val in replace_dict.items():
            file_contents = file_contents.replace(key, str(val))
        file_contents = file_contents.replace("[", "")
        file_contents = file_contents.replace("]", "")
        with open(path_to_inp_file, "w") as file:
            file.write(file_contents)
        # Add oil contents
        with open(oil_file_path, "r") as file:
            oil_contents = file.read()
        with open(path_to_inp_file, "a") as file:
            file.write(oil_contents)
        # Add forecast days dates
        with open(path_to_inp_file, "a") as file:
            n_days = int(simulation["sim_length"] / 24.0) + 1
            file.write(str(n_days) + "\n")
            date_array = np.arange(
                simulation["start_datetime"],
                simulation["start_datetime"] + pd.Timedelta(days=n_days),
                pd.Timedelta(days=1),
            ).astype(datetime.datetime)
            for date in date_array:
                # Convert to 'yyMMdd' format
                formatted_date = date.strftime("%y%m%d")
                file.write(formatted_date)
                file.write("\n")

    def run_medslik_sim(self):
        """
        Run Medslik-II simulation.
        """
        # model directory. Couls be changed, but will remain fixed for the time being.
        model_dir = "src/model/"
        # Compile and start running
        subprocess.run(
            [f"{model_dir}/bin/simulation.exe {self.root_directory}"],
            shell=True,
            check=True,
        )


if __name__ == "__main__":
    # Logging first info
    logger.info("Starting Medslik-II oil spill simulation")
    exec_start_time = datetime.datetime.now()
    logger.info(f"Execution starting time = {exec_start_time}")

    # Config as argument
    parser = ArgumentParser(
        description="Medslik-II # oil spill fate and transport model"
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Path to configuration file",
        required=False,
        default=None,
    )
    args = parser.parse_args()
    config_path = args.config

    if config_path is None:
        config_path = os.path.join("config.toml")
    logger.info("Defining the main object")
    config = Config(config_path).config_dict

    main = MedslikII(config)

    shutil.copy(config_path, os.path.join(main.xp_directory, "config.toml"))

    # Run preprocessing
    logger.info("Starting pre processing ... ")
    preproc = MedslikII.run_preproc(
        main.config,
        main.root_directory,
        main.lon_min,
        main.lon_max,
        main.lat_min,
        main.lat_max,
    )
    logger.info("End of pre processing ...")
    # Create oil file
    oil_path = os.path.join(main.medslik_directory, "data", "oilbase.csv")
    logger.info(f"Reading oil from database {oil_path}")
    oil_value = main.config["simulation"]["oil"][0]
    if isinstance(oil_value, str):
        oil_option = "name"
    else:
        oil_option = "api"
    oil_out_path = read_oilbase(oil_option, oil_value, oil_path, main.xp_directory)
    logger.info(f"Created oil file {oil_out_path}")
    # Create input file
    oilspill_input_file = os.path.join(main.xp_directory, "oilspill.inp")
    main.config_to_input(main.config, oil_out_path, oilspill_input_file)
    logger.info(f"Created input file {oilspill_input_file}")
    # Create parameters file
    parameters_file = os.path.join(main.xp_directory, "parameters.txt")
    if main.config["simulation"]["advanced_parameters"]:
        toml_parameters = main.config["simulation"]["advanced_parameters_path"]
    else:
        toml_parameters = os.path.join(main.medslik_directory, "src", "parameters.toml")
    logger.info(f"Reading parameters file {toml_parameters}")
    if config["simulation"]["advanced_parameters"] == False:
        logger.info("Using custom advanced parameters")
    else:
        logger.info("Using standard parameters")
    main.toml_to_parameters(toml_parameters, parameters_file)
    shutil.copy(toml_parameters, os.path.join(main.xp_directory, "parameters.toml"))
    logger.info(f"Created parameters file {parameters_file}")
    # Run model
    if main.config["run_options"]["run_model"]:
        logger.info("Running Medslik-II simulation")
        main.run_medslik_sim()
    if main.config["run_options"]["postprocessing"]:
        multiple_slick = main.config["simulation"]["multiple_slick"]
        PostProcessing.create_concentration_dataset(
            lon_min=main.lon_min,
            lon_max=main.lon_max,
            lat_min=main.lat_min,
            lat_max=main.lat_max,
            filepath=main.out_directory,
            multiple_slick=multiple_slick,
        )
    if main.config["plot_options"]["plotting"]:
        mplot = MedslikIIPlot(main)
        # Concentration plot
        if _has_pyngl:
            mplot.plot_pyngl(plot_step=1)
        else:
            mplot.plot_matplotlib()
        # mass balance plot
        # mplot.plot_mass_balance()
        try:
            mplot.create_gif()
        except:
            pass
    # Log execution time
    exec_end_time = datetime.datetime.now()
    logger.info(f"Execution ending time = {exec_end_time}")
    logger.info(f"Total execution time = {pd.Timedelta(exec_end_time-exec_start_time)}")
    shutil.move("medslik_run.log", os.path.join(main.root_directory, "medslik_run.log"))
