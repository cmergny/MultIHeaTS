import numpy as np
import matplotlib.pyplot as plt

from multiheats.solvers import ImplicitSolver
from multiheats.old_solvers import OldImplicitSolver, OldValidator


### CLASS


class Validator:
    """Class used to validate the numerical solver solutions.
    It computes analytical solutions for known initial temperatures."""

    def __init__(self, nx) -> None:
        self.nx = nx
        self.temp = np.zeros(nx)
        self.temp[nx // 2 :] = 1
        self.a0 = 1

        self.cond = np.full(self.nx, 0.3)
        self.rho = np.full(self.nx, 1)
        self.cp = np.full(self.nx, 1)

        self.qheat = np.zeros(nx)
        self.eps = 0
        self.L = 1
        self.spaces = np.linspace(0, self.L, self.nx)

        self.alpha = (self.cond / self.rho / self.cp).mean()

    def solve_stepfunction(self, time, n_fourier=50):
        """
        Computes the analytic sol of the heat equa in [0,L] with:
             u(x, 0) = 0 if x<L/2 and 1 if x>L/2
             u'(O) = u'(L) = 0
        Returns: new_temp: array of temperature for t=time, shape = (nx)
        """
        un = np.zeros(n_fourier)
        new_temp = np.zeros(self.nx)

        for ix, x in enumerate(self.spaces):
            # Sum all fourier coefficients
            for n in range(1, n_fourier, 2):
                an = 2 / np.pi / n * (-1) ** ((n - 1) / 2 + 1)
                bn = 0
                wn = n * np.pi * x / self.L  # pulsation
                un[n] = (an * np.cos(wn) + bn * np.sin(wn)) * np.exp(
                    -self.alpha * (n * np.pi / self.L) ** 2 * time
                )
            new_temp[ix] = self.a0 / 2 + un.sum()
        # Properties

        return new_temp


### FUNCTIONS


def compare(spaces, val_temps, sol_temps):
    fig, ax = plt.subplots()
    nt = val_temps.shape[1]
    for it in range(1, nt, nt // 10):
        ax.plot(spaces, val_temps[:, it], color="blue")
        ax.plot(spaces, sol_temps[:, it], color="red")
    ax.set_ylim(0, 1)
    ax.set_xlabel("Depth (m)")
    ax.set_ylabel("Temperature (K)")
    plt.legend(["Analytic", "Solver"])
    plt.show()


def visu(spaces, temps):
    fig, ax = plt.subplots()
    for it in range(temps.shape[1]):
        ax.plot(spaces, temps[:, it], label=f"{it}")
    # ax.set_ylim(0, 1)
    ax.set_xlabel("Depth (m)")
    ax.set_ylabel("Temperature (K)")
    # plt.legend()
    plt.show()


def test_func():
    """
    Test comparaison with analytic solution.
    """
    dt = 1e-4
    nt = 1000
    nx = 100

    val_temps = np.zeros((nx, nt))
    sol_temps = np.zeros((nx, nt))

    valid = Validator(nx)
    solver = ImplicitSolver(valid)
    val_temps[:, 0] = valid.temp
    sol_temps[:, 0] = solver.temp

    # Old solver
    rho = np.ones((nx, nt))
    cp = np.ones((nx, nt))
    k = np.ones((nx, nt)) * 0.3
    old_solver = OldImplicitSolver(
        nx, nt, rho, cp, k, dt, solver.dx[0], valid.L, valid.temp
    )
    old_temps = old_solver.implicit_scheme()

    # Old Validator
    old_valid = OldValidator(nx, nt, rho, cp, k, dt, solver.dx[0], valid.L, valid.temp)
    vold_temps = old_valid.solve_stepfunction()

    # New Solver and Valid
    time = 0
    for it in range(1, nt):
        time += dt
        solver.temp = solver.implicit_scheme(dt)
        sol_temps[:, it] = solver.temp
        val_temps[:, it] = valid.solve_stepfunction(time)

    # visu(valid.spaces, val_temps)
    # visu(valid.spaces, sol_temps)
    # visu(valid.spaces, old_temps)
    compare(valid.spaces, vold_temps, old_temps)


test_func()
