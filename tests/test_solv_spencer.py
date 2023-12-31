import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav
from tqdm import tqdm
import time

import multiheats.constants as cst
from multiheats.solvers import ImplicitSolver


class SpencerModel:
    """docstring for SpencerModel"""

    def __init__(self):
        self.temps = self.import_from_idl("tdarr.dat")
        self.fluxs = self.import_from_idl("flux.dat") / 1e3  # Kt/z
        self.spaces = self.import_from_idl("zarr.dat") / 100  # convert to cm
        self.times = np.linspace(0, 1, self.temps.shape[1])

    def import_from_idl(self, filename):
        """Import array from IDL"""
        file = readsav(f"data/spencer_out/{filename}", verbose=False)
        return file[filename.replace(".dat", "")].transpose()

    def compare(self, v, w, time, filename):
        fig, ax = plt.subplots()
        ax.plot(time, v, ".", color="blue", alpha=0.5, markersize=3)
        ax.plot(time, w, ".", color="red", alpha=0.5, markersize=3)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(filename)
        plt.legend(["Spencer", "Implicit"])
        plt.savefig(filename)


class FakeProfile:
    """
    Fake profile to run the comparaison between Implicit Solver and Spencer Model.
    """

    def __init__(self) -> None:
        self.nx = 101
        self.length = 2.0
        self.spaces = self.regular_space()
        # self.spaces = self.power_spaces(pwr=4)
        self.dx = np.diff(self.spaces)

        self.temp = np.full(self.nx, 89.9130)

        inter = 0.5
        self.rho = bilayer(self.spaces, 800, 2000, inter * 0.5)
        self.cp = bilayer(self.spaces, 600, 1800, inter * 1.0)
        self.inertia = bilayer(self.spaces, 200, 20, inter * 1.5)
        self.cond = self.inertia**2 / self.rho / self.cp

        self.qheat = np.zeros(self.nx)
        self.eps = 1
        self.alb = 0.015

    def regular_space(self):
        spaces = np.linspace(0, self.length, self.nx)
        return spaces

    def power_spaces(self, pwr=2):
        spaces = np.linspace(0, self.length ** (1 / pwr), self.nx) ** pwr
        # plt.plot(spaces, [1 for x in range(self.nx)], ".")
        return spaces


def fake_slr_flux(alb, times, cond):
    """Use a cosinus like function to modelise the solar flux."""
    distance = 9.51 * cst.UA
    period = 79.3 * cst.EARTH_DAY  # day on Japet (in sec)
    solar_flux = (1 - alb) * cst.SOLAR_CST / (distance / cst.UA) ** 2

    surface_bc = -solar_flux * np.cos(2 * np.pi * times / period + np.pi)
    surface_bc[np.where(surface_bc > 0)] = 0
    return surface_bc


def bilayer(x, start, end, inter):
    """Returns an array with a smooth transition
    (tanh) between a start and end values at inter"""
    xprime = (x - inter) / x.max() * 10 * np.pi
    res = (start + end) / 2 + (end - start) / 2 * np.tanh(xprime)
    return res


def plot_error(err, spaces):
    """Unsused"""
    chosen_curves = [3, 100, 423, 34897]
    fig, ax = plt.subplots()
    for it in chosen_curves:
        ax.plot(spaces, err[:, it], label=f"{it}")
    # ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("Depth (m)")
    ax.set_ylabel("Error (K)")
    plt.legend()
    plt.show()


def compare_plots(spencer, temps, prof):
    """Unsused"""
    # chosen_curves = [0, 434, 1288, 3488, -10]
    chosen_curves = [10]
    fig, ax = plt.subplots()
    spencer.spaces[0] = 1e-2
    prof.spaces[0] = 1e-2
    factor = 1e3
    for it in chosen_curves:
        ax.plot(prof.spaces, temps[:, it], ".", color="red")
        # ax.plot(prof.spaces, temps_CN[:, it], ".", color="green")
        ax.plot(spencer.spaces, spencer.temps[:, int(it * factor)], ".", color="blue")
        # ax.plot(prof.spaces[1:], temps[1:, it], ".", color="red")
    # ax.set_xscale("log")
    ax.set_xlabel("Depth (m)")
    ax.set_ylabel("Temperature (K)")
    # ax.set_ylim(30, 150)
    # ax.set_xlim(1e-2, 2)
    plt.show()


def run_solver(name, frame_jump=1):
    """Iterate over time to run solver and return array of temperatures."""

    nday = 5
    nt = int(nday * 10e3 / frame_jump)
    tf = nday * 79.3 * cst.EARTH_DAY
    times = np.linspace(0, tf, nt)

    prof = FakeProfile()
    slr_flux = fake_slr_flux(prof.alb, times, prof.cond)

    nt = times.shape[0]
    temps = np.zeros((prof.nx, nt))
    temps[:, 0] = prof.temp
    dt = np.diff(times)[0]

    start = time.time()

    for it in tqdm(range(1, nt)):
        if name == "Implicit":
            solver = ImplicitSolver(prof)
            prof.temp = solver.implicit_scheme(dt, slr_flux[it])
        # if name == "CN" and it < nt - 1:
        #     solver = CrankNicolson(prof)
        #     prof.temp = solver.CN_scheme(dt, -slr_flux[it], -slr_flux[it + 1])
        temps[:, it] = prof.temp

    duration = (time.time() - start) / nt

    return temps, prof.spaces, duration


def test_spencer(threshold=1.0):
    spencer = SpencerModel()
    solv_temps, _, _ = run_solver("Implicit", frame_jump=1)

    pos = np.arange(0, solv_temps.shape[0])
    pos = np.delete(pos, 1)
    err = abs(spencer.temps - solv_temps[pos, :])

    # plot_error(err, prof.spaces[1:])
    assert err.max() < threshold


if __name__ == "__main__":
    # test_spencer()
    print("Not a main program")
    test_spencer()
