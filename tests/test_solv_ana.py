import numpy as np

# import matplotlib.pyplot as plt

from multiheats.solvers import ImplicitSolver
from multiheats.solvers import CrankNicolson

# import visu_ana as vis


### CLASS


class FakeProfile:
    """
    Fake profile to run the comparaison between Implicit Solver and Analytic solution.
    """

    def __init__(self, nx) -> None:
        self.nx = nx
        self.length = 1
        self.spaces = self.irregular_space()
        self.dx = np.diff(self.spaces)

        self.temp = np.zeros(nx)
        idxs = np.where(self.spaces > self.length / 2)
        self.temp[idxs] = 1

        self.cond = np.full(self.nx, 0.3)
        self.rho = np.full(self.nx, 1)
        self.cp = np.full(self.nx, 1)
        self.qheat = np.zeros(nx)
        self.eps = 0
        self.alpha = (self.cond / self.rho / self.cp).mean()

    def regular_space(self):
        spaces = np.linspace(0, self.length, self.nx)
        return spaces

    def irregular_space(self):
        mean = self.length / 2
        thick = self.length / 5
        spaces = np.linspace(0, self.length, self.nx - 1)
        distances = np.abs(spaces - mean)

        sigmoid = 1 / (1 + np.exp(-(distances - mean) / mean / thick))
        sigmoid /= np.max(sigmoid)

        irr_spaces = np.zeros(self.nx)
        irr_spaces[1:] = np.cumsum(sigmoid) * self.length / np.sum(sigmoid)
        return irr_spaces

    def power_spaces(self):
        pwr = 1.5
        spaces1 = np.linspace(0, (self.length / 2) ** pwr, self.nx // 2) ** (1 / pwr)
        jump = (spaces1[-1] - spaces1[-2]) ** (1 / pwr)
        spaces2 = (
            np.linspace(jump, (self.length / 2) ** (1 / pwr), self.nx // 2) ** pwr
            + self.length / 2
        )
        spaces = np.concatenate((spaces1, spaces2))
        return spaces


### FUNCTIONS


def analytic_step_fct(spaces, time, alpha, n_fourier=50):
    """
    Computes the analytic sol of the heat equa in [0,L] with:
         u(x, 0) = 0 if x<L/2 and 1 if x>L/2
         u'(O) = u'(L) = 0
    Returns: new_temp: array of temperature for t=time, shape = (nx)
    """
    # spaces = np.linspace(0, 1, 100)
    a0 = 1
    length = spaces[-1]
    un = np.zeros(n_fourier)
    new_temp = np.zeros(spaces.shape[0])

    for ix, x in enumerate(spaces):
        # Sum all fourier coefficients
        for n in range(1, n_fourier, 2):
            an = 2 / np.pi / n * (-1) ** ((n - 1) / 2 + 1)
            bn = 0
            wn = n * np.pi * x / length  # pulsation
            un[n] = an * np.cos(wn) + bn * np.sin(wn)
            un[n] *= np.exp(-alpha * (n * np.pi / length) ** 2 * time)
        new_temp[ix] = a0 / 2 + un.sum()
    return new_temp


def test_ana(threshold=0.05, nt=500):
    """
    Runs ImplicitSolver and Validator with same properties.
    Computes the error between the two temperatures over nt iterations.
    Test if mean total error is less than a certain threshold.
    """
    nx = 30

    prof = FakeProfile(nx)
    solver = ImplicitSolver(prof)

    val_temps = np.zeros((nx, nt))
    sol_temps = np.zeros((nx, nt))
    val_temps[:, 0] = prof.temp
    sol_temps[:, 0] = prof.temp

    # Stability Needed for Accuracy
    time = 0
    unstable_factor = 1e3
    dt = unstable_factor * prof.dx.min() ** 2 / (2 * prof.alpha)

    for it in range(1, nt):
        time += dt
        solver.temp = solver.implicit_scheme(dt, 0)

        sol_temps[:, it] = solver.temp
        val_temps[:, it] = analytic_step_fct(prof.spaces, prof.alpha, time, 50)

    err = abs(val_temps - sol_temps)
    assert err.mean() < threshold
    # return err


if __name__ == "__main__":
    err = test_ana()
