# ------------------------------------------------
# MEDSLIK-II oil spill fate and transport model
# ------------------------------------------------
import sys, os
import shutil
import logging
import warnings
import subprocess
import numpy as np
import pandas as pd
from datetime import datetime
from argparse import ArgumentParser
from glob import glob

from src.preprocessing.preprocessing_mdk3 import PreProcessing
from src.postprocessing.postprocessing_mdk3 import create_concentration_dataset
from src.plot.plot_mdk3 import MedslikIIPlot
from src.utils.config import Config
from src.utils.utils import Utils
from download.download_era5_parser import *
from download.download_copernicus_parser import *

import logging

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


class Main:
    """
    This class embeds the MAIN code of medslik-II software.
    """

    def __init__(self, config_path: str) -> None:
        """
        Class constructor given config file path.
        """
        config = Config(config_path).config_dict
        self.config = config
        self.config_path = config_path
        self.root_directory = f"{self.config['input_files']['preproc_path']}{config['simulation']['name']}/"
        self.out_directory = f"{self.root_directory}/out_files/"
        self.out_figures = f"{self.out_directory}figures/"
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

    def apply_aging_effects(self) -> None:
        oilspill = self.config["simulation"]
        age = oilspill["slick_age"]
        if self.n_spill_points > 1:
            oilspill["sim_length"] += age
            oilspill["start_datetime"] -= pd.Timedelta(hours=age)
            oilspill["spill_duration"] = 0
        self.config["simulation"] = oilspill

    def initial_checking(self):

        # checking if the coordinates are on land
        lat = self.config["simulation"]["spill_lat"]
        lon = self.config["simulation"]["spill_lon"]
        coastline_path = self.config["input_files"]["dtm"]["coastline_path"]
        sea = Utils.check_land(lon, lat, coastline_path)

        if sea == 0:
            raise ValueError(
                "Your coordinates lie within land. Please check your values again"
            )

        # checking dates
        dt = Utils.validate_date(self.config["simulation"]["start_datetime"])

        logger.info("No major issues found on dates and oil spill coordimates")

    def data_download_medslik(self):

        config = self.config

        lat_min, lat_max = self.lat_min, self.lat_max
        lon_min, lon_max = self.lon_min, self.lon_max
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
                    f'cp {output_path}*{identifier}*{config["simulation"]["name"]}*.nc {self.root_directory}oce_files/'
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
                    f'cp {output_path}*{identifier}*{config["simulation"]["name"]}*.nc {self.root_directory}met_files/'
                ],
                shell=True,
            )
            subprocess.run(
                [f'rm {output_path}*{identifier}*{config["simulation"]["name"]}*.nc'],
                shell=True,
            )

    def run_preproc(config):

        preproc = PreProcessing(
            config=config,
            input_folder=None,
            preproc_data_path=config["input_files"]["preproc_path"],
        )

        preproc.create_directories()

        # download data if needed
        if config["download"]["download_data"] == True:
            main.data_download_medslik()

        if config["run_options"]["preprocessing"]:

            if config["run_options"]["preprocessing_metoce"]:
                preproc.process_currents()
                preproc.process_winds()

            if config["run_options"]["preprocessing_dtm"]:
                preproc.common_grid()
                preproc.process_bathymetry(
                    config["input_files"]["dtm"]["bathymetry_path"]
                )
                preproc.process_coastline(
                    config["input_files"]["dtm"]["coastline_path"]
                )

            spill_dictionary = {}
            if main.n_spill_points > 1:
                logger.info(
                    f"Starting to write {main.n_spill_points} events of oil spill"
                )
                for i, dur in enumerate(config["spill_rate"]):
                    # obtaining the variables
                    spill_dictionary["simname"] = preproc.simname
                    spill_dictionary["dt_sim"] = main.config["simulation"][
                        "start_datetime"
                    ]
                    spill_dictionary["sim_length"] = preproc.sim_length
                    spill_dictionary["longitude"] = main.config["simulation"][
                        "spill_lon"
                    ][i]
                    spill_dictionary["latitude"] = main.config["simulation"][
                        "spill_lat"
                    ][i]
                    spill_dictionary["spill_duration"] = int(
                        main.config["simulation"]["spill_duration"][i]
                    )
                    spill_dictionary["spill_rate"] = main.config["simulation"][
                        "spill_rate"
                    ][i]
                    spill_dictionary["oil_api"] = main.config["simulation"]["oil"][i]
                    spill_dictionary["number_slick"] = 1
                    preproc.write_config_files(
                        spill_dictionary, separate_slicks=True, s_num=i
                    )

            else:
                logger.info("Writing single slick event")
                spill_dictionary["simname"] = preproc.simname
                spill_dictionary["dt_sim"] = main.config["simulation"]["start_datetime"]
                spill_dictionary["sim_length"] = int(preproc.sim_length)
                spill_dictionary["longitude"] = main.config["simulation"]["spill_lon"][
                    0
                ]
                spill_dictionary["latitude"] = main.config["simulation"]["spill_lat"][0]
                spill_dictionary["spill_duration"] = int(
                    main.config["simulation"]["spill_duration"][0]
                )
                spill_dictionary["spill_rate"] = main.config["simulation"][
                    "spill_rate"
                ][0]
                spill_dictionary["oil_api"] = main.config["simulation"]["oil"][0]
                preproc.write_config_files(spill_dictionary, separate_slicks=False)

            logger.info("Modfying medslik_II.for")
            preproc.process_medslik_memmory_array()
            logger.info("Medslik-II simulation parameters")
            if config["simulation"]["advanced_parameters"] == False:
                logger.info("Using custom advanced parameters")
                preproc.configuration_parameters()
            else:
                logger.info("Using standard parameters")

    def run_medslik_sim(self, simdir, simname, separate_slicks=False):

        # model directory. Couls be changed, but will remain fixed for the time being.
        model_dir = "src/model/"

        day = self.config["simulation"]["start_datetime"].day
        year = self.config["simulation"]["start_datetime"].year
        month = self.config["simulation"]["start_datetime"].month
        hour = self.config["simulation"]["start_datetime"].hour
        minute = self.config["simulation"]["start_datetime"].minute

        output_dir = f"{model_dir}OUT/MDK_SIM_{year}_{month:02d}_{day:02d}_{hour:02d}{minute:02d}_{simname}/."

        # removing old outputes just to be sure
        subprocess.run([f"rm -rf {output_dir}"], shell=True)

        if separate_slicks == False:
            # copy METOCEAN files to MEDSLIK-II installation
            subprocess.run(
                [f"cp {simdir}{simname}/oce_files/*.mrc {model_dir}RUN/TEMP/OCE/"],
                shell=True,
                check=True,
            )
            subprocess.run(
                [f"cp {simdir}{simname}/met_files/*.eri {model_dir}RUN/TEMP/MET/"],
                shell=True,
                check=True,
            )
            # copy bnc files
            subprocess.run(
                [f"cp {simdir}{simname}/bnc_files/* {model_dir}DTM_INP/"],
                shell=True,
                check=True,
            )
            # copy Extract and config files
            subprocess.run(
                [
                    f"cp {simdir}{simname}/xp_files/medslik_II.for {model_dir}RUN/MODEL_SRC/"
                ],
                shell=True,
                check=True,
            )
            subprocess.run(
                [f"cp {simdir}{simname}/xp_files/config2.txt {model_dir}RUN/"],
                shell=True,
                check=True,
            )
            subprocess.run(
                [f"cp {simdir}{simname}/xp_files/config1.txt {model_dir}RUN/"],
                shell=True,
                check=True,
            )
            # Compile and start running
            subprocess.run(
                [f"cd {model_dir}RUN/; sh MODEL_SRC/compile.sh; ./RUN.sh"],
                shell=True,
                check=True,
            )

        else:
            slicks = glob(f"{simdir}{simname}/xp_files/*/")
            for i in range(0, len(slicks)):
                subprocess.run(
                    [f"cp {simdir}{simname}/oce_files/*.mrc {model_dir}RUN/TEMP/OCE/"],
                    shell=True,
                )
                subprocess.run(
                    [f"cp {simdir}{simname}/met_files/*.eri {model_dir}RUN/TEMP/MET/"],
                    shell=True,
                )
                # copy bnc files
                subprocess.run(
                    [f"cp {simdir}{simname}/bnc_files/* {model_dir}DTM_INP/"],
                    shell=True,
                )
                # copy Extract and config files
                subprocess.run(
                    [
                        f"cp {simdir}{simname}/xp_files/medslik_II.for {model_dir}RUN/MODEL_SRC/"
                    ],
                    shell=True,
                )
                subprocess.run(
                    [f"cp {simdir}{simname}/xp_files/config2.txt {model_dir}RUN/"],
                    shell=True,
                )
                subprocess.run(
                    [
                        f"cp {simdir}{simname}/xp_files/slick{i+1}/config1.txt {model_dir}RUN/"
                    ],
                    shell=True,
                )
                # Compile and start running
                subprocess.run(
                    [f"cd {model_dir}RUN/; sh MODEL_SRC/compile.sh; ./RUN.sh"],
                    shell=True,
                    check=True,
                )

        # Send files to case dir and remove temp files
        subprocess.run([f"cp -r {output_dir} {simdir}{simname}/out_files/"], shell=True)
        subprocess.run(
            [f"rm -rf {simdir}{simname}/out_files/MET {simdir}{simname}/out_files/OCE"],
            shell=True,
        )


if __name__ == "__main__":

    logger.info("Starting Medslik-II oil spill simulation")

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
    main = Main(config_path)

    # performing initial checking
    main.initial_checking()

    logger.info("Starting pre processing")
    preproc = Main.run_preproc(main.config)
    logger.info("End of pre processing")

    if main.config["run_options"]["run_model"]:
        logger.info("Running Medslik-II simulation")
        main.run_medslik_sim(
            "cases/", main.config["simulation"]["name"], separate_slicks=False
        )

    if main.config["run_options"]["postprocessing"]:
        multiple_slick = main.config["simulation"]["multiple_slick"]
        create_concentration_dataset(
            lon_min=main.lon_min,
            lon_max=main.lon_max,
            lat_min=main.lat_min,
            lat_max=main.lat_max,
            filepath=main.out_directory,
            multiple_slick=multiple_slick,
        )

    if main.config["plot_options"]["plotting"]:
        mplot = MedslikIIPlot(main)
        available_plots = ["pyNGL", "Matplotlib"]
        plot_product = available_plots[0]
        match plot_product:
            case "pyNGL":
                mplot.plot_pyngl(plot_step=12, concentration_range=None)
            case "Matplotlib":
                mplot.plot_matplotlib()
            case _:
                print("Plotting using default pyNGL")
                mplot.plot_pyngl(plot_step=12, concentration_range=None)
        try:
            mplot.plot_mass_balance()
        except:
            pass
        mplot.create_gif()

    shutil.copy("medslik_run.log", f"{main.out_directory}medslik_run.log")

    if config_path is None:
        shutil.copy("config.toml", f"{main.out_directory}config.toml")
