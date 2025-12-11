import numpy as np
import astropy
from datetime import datetime, timedelta
from astroquery.jplhorizons import Horizons
import warnings
from astropy.units import UnitsWarning

import multiheats.constants as cst


def get_ephemerides(target=None, start=None, stop=None, dt=None, ephfilepath=None, write=True):
    """
    From an observer, start, stop dates,
    Returns the ephemerides of the given target at each timestep.
    Check if
    """
    warnings.filterwarnings("ignore", category=UnitsWarning)
    # Either request SPICE for eph data
    if ephfilepath is None:
        print("Requesting JPL Horizons for ephemerides data...")
        timestep = f"{int(dt/60)}min"  # Horizons format
        targetid = cst.NAIF_DIC[targetname]
        obj = Horizons(
            id=targetid,
            location="@sun",  # Location refers to location of the observer
            epochs={"start": start, "stop": stop, "step": timestep},
        )
        eph = obj.ephemerides()

        if write:
            filename = f"{target}_{start}_{stop}"
            eph.write(f"../../ephemerides/{filename}.fits", format="fits")
            print(f"Wrote SPICE data in {filename}.")

    # Or import from provided file
    else:
        print(f"Using ephemerides from file {ephfilepath}")
        eph = astropy.table.Table.read(ephfilepath)
    return eph


def retrieve_slr_fluxs(alb, lat, lon, eph):
    """Get the solar flux array for a given point, based on ephemerides data."""
    print(f"Retrieving solar flux from {eph['targetname'][0]} at lat {lat} and long {lon}.")
    # Retrieve parameters
    lons_sub = np.array(eph["PDSunLon"], dtype=float)  # Subsolar point longitude
    distances = np.array(eph["r"], dtype=float) * cst.UA  # Distance to the sun
    vis = eph["sat_vis"]  # Tells if primary body occults the target!
    nt = distances.size  # nbr of iterations

    # Compute Solar Flux
    slr_fluxs = np.zeros(nt)
    for it in range(nt):
        if vis[it] in ["U", "P"]:  # If during an eclipse then no flux
            slr_fluxs[it] = 0  # /!\ chose that partial "P" = 0
        else:
            slr_fluxs[it] = compute_solar_flux(alb, lat, lon, lons_sub[it], distances[it])
    return slr_fluxs


def compute_solar_flux(alb, lat, lon, lon_sub, distance):
    """
    From the solar incidence angle, distance and Albedo
    Returns solar flux
    """
    # Compute vector normals
    rsurf = get_vector(np.radians(lat), np.radians(lon))
    rsun = get_vector(0, np.radians(lon_sub))
    # Compute Incidence angle
    cos_inc = np.dot(rsun, rsurf)
    cos_inc = np.clip(cos_inc, 0, 1)  # Night side → 0
    # Solar Flux
    return (1 - alb) * cst.SOLAR_CST / (distance / cst.UA) ** 2 * cos_inc


def get_vector(lat, long):
    """Return vector in cartesian coordinates form lat and long."""
    xsurf = np.cos(lat) * np.cos(long)
    ysurf = np.cos(lat) * np.sin(long)
    zsurf = np.sin(lat)
    return np.array([xsurf, ysurf, zsurf])
    
def find_idx_period(times_utc, heliorb_per):
    """Find index when times becomes periodic"""
    target_time = times_utc[0] + timedelta(seconds=heliorb_per)
    return np.where(times_utc > target_time)[0][0]
        
def find_equilibrium(prof, solver, norbs, periodic_fluxs, dt, heliorb_per, skip_eq):
    """Compute multiple orbits to find equilibrium temperature"""
    if not skip_eq:
        # Initalize
        ntper = periodic_fluxs.size
        temps4eq = np.zeros((ntper, prof.nx))
        
        for iorb in range(norbs): # Each orbit
            print(f"Searching for equilibrium {iorb+1}/{norbs} orbit...")
            for it in range(ntper): # All iter of one orbit
                solver.temp = solver.implicit_scheme(dt, periodic_fluxs[it])
                solver.need_update = False
                
                if iorb == norbs-1: # Save last orbit
                    temps4eq[it, :] = solver.temp  
                    
        # Take Mean Temperature over relevant depths
        bot = 3 * prof.thermal_skin(heliorb_per)[0]
        ibot = np.argmin(np.abs(prof.spaces - bot))
        temp_eq = np.mean(temps4eq[:, :ibot]) # Up to bottom limit
    else:
        print("Skipped equibrilium search")
    return temp_eq


if __name__ == "__main__":

    target = "Ganymede"  # See NAIF ID
    start = "2012-11-13"  # UTC format YYYY-MM-DD
    stop = "2025-11-13"  # UTC format YYYY-MM-DD
    dt = 2 * 3600  # (s) Timestep
    # Point
    alb = 0.5
    lat = 10
    lon = 20
    ephfilepath = "../../ephemerides/eph_503_2012-11-13_2025-11-13.fits"
    # ephfilepath = None

    eph = get_ephemerides(target, start, stop, dt, ephfilepath)  # Ephemerides
    slr_fluxs = retrieve_slr_fluxs(alb, lat, lon, eph)

    import matplotlib.pyplot as plt

    # times = eph["datetime_str"]
    nt = slr_fluxs.size
    times = np.arange(nt)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(times, slr_fluxs, "-o", markersize=2)
    ax.set_xlabel("Time (UTC)")
    ax.set_ylabel("Solar Flux (SI)")
    plt.tight_layout()
    # ax.grid(True)
    plt.show()
