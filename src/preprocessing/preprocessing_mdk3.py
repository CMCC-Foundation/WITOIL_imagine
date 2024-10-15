import sys, os
import numpy as np
import pandas as pd
import xarray as xr
import logging
import shutil
import geopandas as gpd
import subprocess
from glob import glob
from numpy.typing import NDArray
from datetime import datetime, timedelta, date

from src.utils.utils import Utils

import logging

logger = logging.getLogger(__name__)


class PreProcessing:
    """
    This class embeds functions for pre-processing, with a user-interface purpose.
    """

    def __init__(self, config: dict, exp_folder: str, domain=list[float]) -> None:
        """
        Class constructor.
        """
        self.exp_folder = exp_folder
        os.makedirs(exp_folder, exist_ok=True)
        # - start_date and end_date: start and end data of the pre-processing
        self.config = config
        self.day = config["simulation"]["start_datetime"].day
        self.year = config["simulation"]["start_datetime"].year
        self.month = config["simulation"]["start_datetime"].month
        self.hour = config["simulation"]["start_datetime"].hour
        self.minute = config["simulation"]["start_datetime"].minute
        self.sim_length = config["simulation"]["sim_length"]
        start_date = date(self.year, self.month, self.day)
        sim_length_days = int(self.sim_length / 24.0) + 1
        end_date = start_date + timedelta(days=sim_length_days)
        self.start_date = start_date
        self.end_date = end_date
        self.date_list = pd.date_range(start_date, end_date, freq="D")
        self.grid = None
        self.domain = domain

    def create_directories(self):
        os.makedirs(f"{self.exp_folder}/oce_files", exist_ok=True)
        os.makedirs(f"{self.exp_folder}/met_files", exist_ok=True)
        os.makedirs(f"{self.exp_folder}/bnc_files", exist_ok=True)
        os.makedirs(f"{self.exp_folder}/detections", exist_ok=True)

    def process_currents(self, oce_path: str = None):

        logger.info("Pre processing currents")
        lon_min, lon_max, lat_min, lat_max = self.domain
        # opening all files in the directory and concatenating them automatically through open_mfdataset
        if oce_path is None:
            oce_path = f"{self.exp_folder}/oce_files/"
        oce_path = os.path.join(oce_path, "*.nc")
        if glob(oce_path) == []:
            oce_path = f"{self.exp_folder}/oce_files/*.nc"
        concat = xr.open_mfdataset(oce_path, combine="nested", engine="netcdf4")
        concat = concat.drop_duplicates(dim="time", keep="last")
        if self.config["input_files"]["set_domain"]:
            concat = concat.sel(
                lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max)
            )
        # Interpolating the values in time, transforming it from daily to hourly values
        concat = concat.resample(time="1h").interpolate("linear")
        Utils.write_mrc(concat, exp_folder=self.exp_folder)
        self._copy_nc_files(oce_path, f"{self.exp_folder}/oce_files/")

    def process_winds(self, met_path: str = None):

        logger.info("Pre processing winds")
        lon_min, lon_max, lat_min, lat_max = self.domain
        # opening all files in the directory and concatenating them automatically through open_mfdataset
        if met_path is None:
            met_path = f"{self.exp_folder}/met_files/*.nc"
        met_path = os.path.join(met_path, "*.nc")
        if glob(met_path) == []:
            met_path = f"{self.exp_folder}/met_files/*.nc"
        concat = xr.open_mfdataset(met_path, combine="nested", engine="netcdf4")
        concat = concat.drop_duplicates(dim="time", keep="first")
        # Interpolating the values in time, transforming it from daily to hourly values
        concat = concat.resample(time="1h").interpolate("linear")
        # Handle longitude and latitude
        concat["lon"] = xr.where(
            concat["lon"] > 180, concat["lon"] - 360, concat["lon"]
        )
        concat = concat.sortby("lat")
        concat = concat.sortby("lon")
        # Iterating at each hour to generate the .eri files
        for date in self.date_list:
            # Call write eri function located in medslik.utils file
            Utils.write_eri(concat, date, exp_folder=self.exp_folder)
        self._copy_nc_files(met_path, f"{self.exp_folder}/met_files/")

    def _copy_nc_files(self, src_files: str, dst_dir: str) -> None:
        # Use glob to find all .nc files in the source directory
        nc_files = glob(src_files)
        # Loop through each file
        for src_file in nc_files:
            # Get the destination file path
            dst_file = os.path.join(dst_dir, os.path.basename(src_file))
            # Check if the destination file already exists
            if not os.path.exists(dst_file):
                # Copy the file if it doesn't exist
                shutil.copy(src_file, dst_file)

    def process_bathymetry(self, gebco):

        grid = self.grid
        gebco = xr.open_dataset(gebco)

        try:
            grid = grid.rename({"nav_lat": "lat", "nav_lon": "lon"})
        except:
            pass

        # interpolation on medslik grid
        med = gebco.interp(lon=grid.lon.values.tolist(), lat=grid.lat.values.tolist())
        # converting from begative depths to positive
        med["elevation"] = med.elevation * -1
        # filling nan to -9999 as expected by medslik
        med = med.fillna(9999)

        # Convert bathymetry to MDK-II
        mdk_z = []

        llat, llon = len(med.lat), len(med.lon)

        for i in reversed(range(0, llat)):
            for j in range(0, llon):
                rec = med.isel(lat=i, lon=j)
                mdk_z.append(rec.elevation.values.max())

        mdk_z = np.array(mdk_z)
        land_mask = np.where(mdk_z <= 0)
        mdk_z[land_mask] = 9999

        BathFile = open(f"{self.exp_folder}/bnc_files/dtm.bath", "w")
        BathFile.write(
            "MEDSLIK-II compatible bathymetry file. Degraded resolution based on GEBCO 30''\n"
        )
        BathFile.write(
            "%-7.4f %-7.4f %-7.4f %-7.4f \n"
            % (
                np.min(grid.lon.values),
                np.max(grid.lon.values),
                np.min(grid.lat.values),
                np.max(grid.lat.values),
            )
        )
        BathFile.write("%d %d \n" % (llon, llat))
        np.savetxt(BathFile, mdk_z, fmt="%04.0f")
        BathFile.close()

    def process_coastline(self, gshhs):

        grid = self.grid

        try:
            grid = grid.rename({"nav_lat": "lat", "nav_lon": "lon"})
        except:
            pass

        # 1 degree buffer to collect more coastline
        buffer = 1
        # defining minimum and maximum of coordinates box search
        xmin = grid.lon.min() - buffer
        xmax = grid.lon.max() + buffer
        ymin = grid.lat.min() - buffer
        ymax = grid.lat.max() + buffer

        shp = gpd.read_file(gshhs)

        # Cropping to a smaller area
        shp = shp.cx[xmin:xmax, ymin:ymax]

        # shp with linestring instead of polygons
        shp["geometry"] = shp.geometry.boundary

        # Cropping the selected linestrings to the same bounding box
        shp = shp.clip_by_rect(
            xmin.values.max(), ymin.values.max(), xmax.values.max(), ymax.values.max()
        )

        # Removing empty geometries
        shp = shp[~shp.is_empty]

        # Transforming it back again to geodataframe
        shp = gpd.GeoDataFrame(geometry=shp)

        # removing any multiline strings left on the shapefile
        shp = shp.explode(index_parts=True)

        # writes the first line of the .map file. It should contain the # of "isles"
        CoastFile = open(f"{self.exp_folder}/bnc_files/dtm.map", "w")
        CoastFile.write("%-4.0f \n" % (len(shp)))
        iTargetSites = []
        iTargetSite = np.array([0.0, 0.0, 0.0, 0.0])

        for i, polygon in shp.iterrows():

            # extracts the island
            pol = polygon.geometry

            # Extract the exterior coordinates of the polygon
            # exterior_coords = list(pol.exterior.coords)
            exterior_coords = list(pol.coords)

            # prints the length of the island
            CoastFile.write("%-4.0f %1.0f \n" % (len(exterior_coords), 0))
            # prints data related to the island
            for segment in exterior_coords:
                CoastFile.write("%-8.5f %-6.5f \n" % (segment[0], segment[1]))

    def common_grid(self):
        """
        function to crop files based on this presented common grid.

        Default uses current grid for the area of interest
        """
        lon_min, lon_max, lat_min, lat_max = self.domain
        grid = glob(f"{self.exp_folder}/oce_files/*.nc")[0]
        grid = xr.open_dataset(grid)
        if self.config["input_files"]["set_domain"]:
            grid = grid.sel(lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max))
        self.grid = grid

    def write_initial_input(self, coordinates):
        """
        This script considers the values passes on config.toml on area vertex values.

        Please be careful on writing custom domain. Take into considerations values ans the syntax described.

        1st list - Comprehends all slicks
        2nd list - Comprehends all vertex for a single slick
        3rd list - Latitude and Longitude coordinate for a single vertex within one slick

        Writing initial slick countour to be used inside simulation.for
        """

        len_data_points = []
        data_points = ""

        for slick in coordinates:
            # Generate the formatted coordinates
            formatted_coordinates = "\n".join(
                [f"{lat:.5f} {lon:.5f}" for lat, lon in slick]
            )

            # Add the first point as final to close polygon
            closing_point = f"{slick[0][0]:.5f} {slick[0][1]:.5f}\n"  # Replace with the actual value if needed

            # write the coordinates in data_points variable
            data_points = data_points + formatted_coordinates + "\n" + closing_point

            # save the length of data points
            len_data_points.append(len(slick) + 1)

        # Prepare the header
        date = f"{self.year}/{self.month:02d}/{self.day:02d} {self.hour:02d}:{self.minute:02d}"
        num_data_points = sum(
            len_data_points
        )  # sum of all slicks, considering the closing point
        header = (
            f"{date}\n{num_data_points:>2}        Number of data points\n  lat     lon"
        )

        # Full output combining the header, coordinates, and the final point
        output = f"{header}\n{data_points}"

        # Specify the file path where you want to save the file
        file_path = f"{self.exp_folder}/xp_files/contour_slick.csv"

        # Write the output to a text file
        with open(file_path, "w") as file:
            file.write(output)

    def process_initial_shapefile(self):
        """

        This script will generate the initial conditions of an area spill for Medslik-II starting from a shapefile.

        If the shapefile has one os more slicks, the script will consider and write them.

        """

        coastline = gpd.read_file(self.config["input_files"]["dtm"]["coastline_path"])

        #calling the domain
        lon_min, lon_max, lat_min, lat_max = self.domain

        #cropping the coastline polygon in the area of interest
        coastline = coastline.clip_by_rect(xmin=lon_min,xmax=lon_max,ymin=lat_min,ymax=lat_max)

        #small buffer on coastline to ensure no oil inside coast
        b_coastline = coastline.buffer(0.0035)

        obs = gpd.read_file(self.config["input_files"]["shapefile"]["shape_path"])

        ######## OBTAINING THE SLICKS FOR ALL GEOMETRIES IN THE SHAPEFILE ########

        # Explodes different polygons in different entities, so multiple slicks can be generated
        # shape_obs = obs.explode().reset_index()

        # Using convex hull to simplufy the shapefile and start simulation from it
        shape_obs = obs.convex_hull

        #removing possible oil spill areas inside coastline
        crop = gpd.GeoDataFrame(geometry=shape_obs.geometry.difference(b_coastline.geometry))

        crop = crop.dropna()

        # get the list of coordinates
        coords = crop.get_coordinates()

        coords = coords[["y", "x"]]

        coords.columns = ["lat", "lon"]

        output = f"{self.year}/{self.month:02d}/{self.day:02d} {self.hour:02d}:{self.minute:02d}\n"
        output += f"{len(coords)}        Number of data points\n"
        output += "  lat     lon\n"

        # Iterate over the rows to format each lat/lon to 5 decimal places
        for _, row in coords.iterrows():
            output += f"{row['lat']:.5f} {row['lon']:.5f}\n"

        file_path = f"{self.exp_folder}/xp_files/contour_slick.csv"

        # Write the output to a text file
        with open(file_path, "w") as file:
            file.write(output)

if __name__ == "__main__":
    pass
