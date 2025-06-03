### IMPORTS

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

import multiheats as mheats
import multiheats.visualise as vis
import multiheats.constants as cst


# PARAMETERS

nx = 100  # Grid points
xmin, xmax = 0, 2  # depth limits (m)
alb = 0.2  # Albedo
eps = 1.0  # Emissivity
nday = 200  # Nbr of days
step_per_day = int(1e2)  # Points per day
distance = 9.51 * cst.UA  # Distance to sun (m)
period = 79.3 * cst.EARTH_DAY  # Diurnal period (s)

# TOP
cond_top = 0.01
rho_top = 917.0
cp_top = 839
# BOTTOM
cond_bot = cond_top / 2
rho_bot = rho_top / 2
cp_bot = cp_top / 2

### MAIN

if __name__ == "__main__":

    prof = mheats.create_profiles.Profile(nx, eps, xmin, xmax, power=3)
    # Ini with cond_top just to get interf
    prof.cond = cond_top
    prof.rho = rho_top
    prof.cp = cp_top
    # Interface
    thermal_skin = prof.thermal_skin(period)
    interface = 2 * thermal_skin  # (m)
    # Create Bilayer Profile
    prof.bilayer_prof(cond_top, rho_top, cp_top, cond_bot, rho_bot, cp_bot, interface)

    # Compute time variables
    nt = nday * step_per_day
    tf = nday * period
    times = np.linspace(0, tf, nt)
    dt = np.diff(times).mean()

    # Guess Eq. Temp
    pertimes = np.arange(0, period, dt)  # time over the full cycle
    temp_eq = mheats.create_fluxs.get_eq_temp(pertimes, distance, period, alb, eps)
    # Set Initial Temp
    prof.temp = np.full(prof.nx, temp_eq)
    temps = np.zeros((nt, prof.nx))

    # Run Solver
    solver = mheats.solvers.ImplicitSolver(prof)
    print("Computing temperature evolution...")
    for it in tqdm(range(times.size)):
        temps[it] = prof.temp
        slr_flux = mheats.create_fluxs.fake_slr_flux(alb, times[it], distance, period)
        solver.temp = solver.implicit_scheme(dt, slr_flux)
        prof.temp = solver.temp
        solver.need_update = False

    print("Visualisation")
    # vis.plot_temp(prof.spaces, temps, it=10, interf=prof.interf)
    vis.plot_multi_temp(prof.spaces, temps, n_curves=7)
    plt.savefig("temp_profiles.png")
    anim = vis.animate_function(
        prof.spaces, temps, interf=prof.interf, step=3, frames=400, save=False
    )
    plt.show()
