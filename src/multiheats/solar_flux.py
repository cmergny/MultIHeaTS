import pandas as pd
import numpy as np
from scipy.io import readsav

# import matplotlib.pyplot as plt

import multiheats.constants as cst


class SurfFlux:
    """
    Import solar flux and albedos and generature boundary condition.
    Right now works for Japet files only, but will be expanded in the future.

    """

    def __init__(self) -> None:
        albedo_path = "../../data/Blackburn_Albedo_Iapetus.xls"
        flux_path = "../../data/Absorbed_flux_SATURN_IAPETUS.sav"
        self.import_flux(flux_path)
        self.import_albedos(albedo_path)

    def import_flux(self, filepath):
        """
        Read flux on Japet from idl saved.
        Return array of fluxs for each lat/long.
        """
        file = readsav(filepath)
        self.total_fluxs = file["tab_flux_sol"]
        self.flux_longs = file["tab_long"]
        self.flux_lats = file["tab_lat"]
        self.times = file["tab_time"]

    def import_albedos(self, filepath):
        """
        Read xls file containing albdos on Japet
        Return array of albedos, and vectors containg lats and longs.
        """
        df = pd.read_excel(filepath, header=0)
        albedos = df.to_numpy()
        # First line is longs
        longs = list(df)
        longs[0] = 0
        self.alb_longs = np.array(longs)
        # First column is lats
        self.alb_lats = albedos[:, 0]
        self.albedos = albedos[:, 1:]

    def get_local_fluxs(self, lat, long):
        ilat = find_nearest(self.flux_lats, lat)
        ilong = find_nearest(self.flux_longs, long)
        return self.total_fluxs[ilat, ilong]

    def get_local_albedos(self, lat, long):
        ilat = find_nearest(self.alb_lats, lat)
        ilong = find_nearest(self.alb_longs, long)
        return self.albedos[ilat, ilong]

    def get_eq_temp(self, lat, long, eps):
        """Get close to eq temp"""
        alb = self.get_local_albedos(lat, long)

        mean_slr_flux = self.get_local_fluxs(lat, long).mean()
        mean_surf_flux = (1 - alb) * mean_slr_flux

        temp_eq = (mean_surf_flux / eps / cst.SIGMA) ** (1 / 4)
        return temp_eq

    def get_solar_fluxs(self, lat, long):
        """
        Get flux array for a given lat and long.
        in W/m2
        """
        flux = self.get_local_fluxs(lat, long)
        self.alb = self.get_local_albedos(lat, long)

        return (1 - self.alb) * flux


def find_nearest(array, value, verbose=False):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if verbose:
        print(f"Input value {value}, taken value {array[idx]}.")
    return idx


if __name__ == "__main__":

    print("This program is used for imports only.")
    # lat, long = 13, 78

    # surf = SurfFlux()
    # times = surf.times
    # fluxs = np.zeros(times.shape[0])

    # for it, time in enumerate(times):
    #     slr_flux, alb = surf.get_flux(time, lat, long)
    #     fluxs[it] = slr_flux

    # fig, ax = plt.subplots()
    # ax.plot(times, fluxs, ".")
    # ax.set_xlabel("Time")
    # ax.set_ylabel("Flux")
    # plt.title(f"Lat {lat}, Long {long}")
    # plt.show()
