### IMPORTS

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

from multiheats.solvers import ImplicitSolver
from multiheats.create_profiles import Profile
import multiheats.visualise as vis
import multiheats.constants as cst
import multiheats.create_fluxs as create_fluxs


# PARAMETERS

nx = 100  # Grid points
xmin, xmax = 0, 1  # depth limits (m)
alb = 0.2  # Albedo
eps = 1.0  # Emissivity
nday = 200  # Nbr of days
step_per_day = int(1e2)  # Points per day
distance = 9.51 * cst.UA  # Distance to sun (m)
period = 79.3 * cst.EARTH_DAY  # Diurnal period (s)

prof = Profile(nx, eps, xmin, xmax, power=3)
# TOP
cond_top = 0.01
rho_top = 917.0
cp_top = 8390
# BOTTOM
cond_bot = cond_top / 2
rho_bot = rho_top / 2
cp_bot = cp_top / 2
# Interface
thermal_skin = prof.thermal_skin(cond_top, rho_top, cp_top, period)
interface = 2 * thermal_skin  # (m)


### MAIN

if __name__ == "__main__":
    # Create Bilayer Profile
    prof.bilayer_prof(cond_top, rho_top, cp_top, cond_bot, rho_bot, cp_bot, interface)

    # Compute time variables
    nt = nday * step_per_day
    tf = nday * period
    times = np.linspace(0, tf, nt)
    dt = np.diff(times).mean()

    # Create Flux Data
    solar_fluxs = create_fluxs.fake_slr_flux(alb, times, distance, period)
    # Guess Initial Temp
    temp_eq = create_fluxs.get_eq_temp(-solar_fluxs, alb, eps)
    prof.temp = np.full(prof.nx, temp_eq)
    temps = np.zeros((nt, prof.nx))

    # Run Solver
    solver = ImplicitSolver(prof)
    print("Computing temperature evolution...")
    for it in tqdm(range(times.size)):
        temps[it] = prof.temp
        prof.temp = solver.implicit_scheme(dt, solar_fluxs[it])
        solver.temp = prof.temp
        solver.need_update = False

    print("Visualisation")
    # vis.plot_temp(prof.spaces, temps, it=10, interf=prof.interf)
    vis.plot_multi_temp(prof.spaces, temps, n_curves=7)
    plt.savefig("temp_profiles.png")
    anim = vis.animate_function(
        prof.spaces, temps, interf=prof.interf, step=3, frames=400, save=False
    )
    plt.show()
