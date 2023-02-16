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
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.animation import FuncAnimation

from solvers import ImplicitSolver
from create_profiles import Profile
from solar_flux import SurfFlux
import visualise as vis

### MAIN


def animate_function(spaces, temps, interf=0, step=1, frames=None, save=False):
    """Plot an animation of the temperature with time"""
    frames = temps.shape[0] // step if frames is None else frames
    fig, ax = plt.subplots()
    (line,) = ax.plot([])
    ax.axvline(x=interf, alpha=0.3, linestyle="-", color="black")

    def animate(it):
        line.set_data((spaces, temps[it * step]))
        # ax.set_title(f"{it}")
        return line

    ax.set_ylabel("Temperature (K)")
    ax.set_xlabel("Depth (m)")
    ax.set_xlim(spaces.min(), spaces.max())
    ax.set_ylim(temps.min(), temps.max())
    ax.set_xscale("symlog")
    anim = FuncAnimation(fig, animate, frames=frames, interval=1, repeat=False)
    plt.title("Temperature Evolution")
    return anim


if __name__ == "__main__":

    print("Creating surface profile...")
    prof = Profile()
    prof.bilayer_prof()
    # prof.monolayer_prof()

    # Get Initial Temperature
    surf = SurfFlux()
    temp_eq = surf.get_eq_temp(prof.lat, prof.long, prof.eps)
    prof.temp = np.full(prof.nx, temp_eq)

    times = surf.times[:]
    nt = times.shape[0]
    dt = np.diff(times)[0]
    temps = np.zeros((nt, prof.nx))

    print("Computing temperature evolution...")
    for it in tqdm(range(times.shape[0])):
        solver = ImplicitSolver(prof)
        solver.solar_flux = -surf.get_flux(times[it], prof.lat, prof.long)
        prof.temp = solver.implicit_scheme(dt)

        temps[it] = prof.temp

    print("Visualisation")
    # Visualisation
    it = 300
    # vis.plot_temp(prof.spaces, temps, it)
    # vis.plot_multi_temp(prof.spaces, temps, n_curves=10)

    anim = animate_function(prof.spaces, temps, interf=prof.interf, step=5, frames=350)
    plt.show()

    tex_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 13,
        "font.size": 11,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
    }
    plt.rcParams.update(tex_fonts)

    progress_callback = lambda i, n: print(f"Saving frame {i} of {n}")
    anim.save("animation.gif", progress_callback=progress_callback)
