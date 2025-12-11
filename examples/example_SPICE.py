### IMPORTS

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import datetime

import multiheats as mheats
import multiheats.visualise as vis
import multiheats.constants as cst


# PARAMETERS

# Solar flux
target = "503"  # See NAIF ID
start = "2012-11-13"  # UTC format YYYY-MM-DD
stop = "2025-11-13"  # UTC format YYYY-MM-DD
dt = 2 * 3600  # (s) Timestep
heliorb_per = 11.8618 * cst.YEAR  # Heliocentric orbit

# Point
alb = 0.5
lat = 10
lon = 20
ephfilepath = "../ephemerides/eph_503_2012-11-13_2025-11-13.fits"
# ephfilepath = None

# SPACE
nx = 100  # Grid points
xmin, xmax = 0, 50  # depth limits (m)
xpower = 5

# TOP
cond = 0.01
rho = 917.0
cp = 839
eps = 0.9  # Emissivity


def find_idx_period(times_utc, heliorb_per):
    """Find index when times becomes periodic"""
    target_time = times_utc[0] + datetime.timedelta(seconds=heliorb_per)
    return np.where(times_utc > target_time)[0][0]


### MAIN

if __name__ == "__main__":

    # Create Monolayer Material
    prof = mheats.create_profiles.Profile(nx, eps, xmin, xmax, power=xpower)
    prof.cond = cond
    prof.rho = rho
    prof.cp = cp
    prof.monolayer_prof(cond, rho, cp)

    # Retrieve solar flux
    eph = mheats.query_spice.get_ephemerides(target, start, stop, dt, ephfilepath)  # Ephemerides
    slr_fluxs = mheats.query_spice.retrieve_slr_fluxs(alb, lat, lon, eph)
    distances = eph["r"]
    nt = slr_fluxs.size

    # Initialize Temperature
    prof.temp = (np.mean(slr_fluxs) / eps / cst.SIGMA) ** (1 / 4)  # Equilibrium temp
    temps = np.zeros((nt, nx))
    print(f"Initial Temperature set to {prof.temp:.1f} K.")

    # Run Solver
    solver = mheats.solvers.ImplicitSolver(prof)
    solver.temp = np.full(nx, prof.temp)
    print("Computing temperature evolution...")

    # Run multiple orbit to find Equilibrium
    norbs = 20
    times_utc = np.array(
        [datetime.datetime.strptime(t, "%Y-%b-%d %H:%M") for t in eph["datetime_str"]]
    )
    itper = find_idx_period(times_utc, heliorb_per)
    for iorb in range(norbs):
        print(f"Searching for equilibrium {iorb+1}/{norbs} orbit...")
        for it in range(itper):
            solver.temp = solver.implicit_scheme(dt, slr_fluxs[it])
            solver.need_update = False
    solver.temp = np.full(nx, solver.temp[-1])

    print("Starting evolution")
    for it in tqdm(range(nt)):
        solver.temp = solver.implicit_scheme(dt, slr_fluxs[it])
        solver.need_update = False
        temps[it] = solver.temp

    print("Visualisation")
    # vis.plot_temp(prof.spaces, temps, it=10, interf=prof.interf)
    vis.plot_multi_temp(prof.spaces, temps, n_curves=7, log=True)
    plt.savefig("temp_profiles.png")
    anim = vis.animate_function(
        prof.spaces, temps, interf=prof.interf, step=3, frames=400, save=False
    )
    plt.show()
