#!/usr/bin/env python

"""
Implicit solver of the heat equation.

From a given set of heat properties and boundary conditions, 
solves the heat equation using an implicit scheme.
Surface condition is computed using the incoming solar flux 
and the black body emission.
Output an array containing a temperature profile for each time step.
Mathematical derivation of the solver can be found at [contact me by email].
"""

__author__ = "Cyril Mergny"
__credits__ = ["Cyril Mergny"]
__license__ = "Gnu GPL"
__email__ = "cyril.mergny@universite-paris-saclay.fr"

### IMPORTS

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

from multiheats.solvers import ImplicitSolver
from multiheats.create_profiles import Profile
from multiheats.solar_flux import SurfFlux
import multiheats.visualise as vis

### MAIN


if __name__ == "__main__":

    print("Creating surface profile...")
    prof = Profile()
    # prof.bilayer_prof()
    prof.monolayer_prof()

    # Get Initial Temperature
    surf = SurfFlux()
    temp_eq = surf.get_eq_temp(prof.lat, prof.long, prof.eps)
    prof.temp = np.full(prof.nx, temp_eq)

    times = surf.times[:]
    nt = times.shape[0]
    dts = np.diff(times)
    temps = np.zeros((nt, prof.nx))

    solar_fluxs = -surf.get_solar_fluxs(prof.lat, prof.long)

    solver = ImplicitSolver(prof)
    print("Computing temperature evolution...")
    for i in range(25):
        for it in tqdm(range(times.shape[0] - 1)):
            prof.temp = solver.implicit_scheme(dts[it], solar_fluxs[it])
            temps[it] = prof.temp
            solver.temp = prof.temp
            solver.need_update = False

    print("Visualisation")
    it = 10
    # vis.use_latex()
    # vis.plot_temp(prof.spaces, temps, it, interf=prof.interf)
    # vis.plot_multi_temp(prof.spaces, temps, n_curves=10)
    anim = vis.animate_function(
        prof.spaces, temps, interf=prof.interf, step=5, frames=400, save=False
    )
    # plt.savefig("hey.png")
    plt.show()
