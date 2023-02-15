"""
Implicit solver of the heat equation.
"""

### IMPORTS
import numpy as np

import constants as cst

### CLASS


class ImplicitSolver:
    """
    Solver class which uses an implicit
    method to solve the heat equation.
    """

    def __init__(self, prof) -> None:
        # Parameters
        self.name = "implicit"
        # Properties
        self.temp = prof.temp
        self.cond = prof.cond
        self.rho = prof.rho
        self.cp = prof.cp
        self.qheat = prof.qheat
        self.eps = prof.eps
        self.nx = prof.temp.shape[0]
        self.solar_flux = 0
        self.dx = np.diff(prof.spaces)

    def implicit_scheme(self, dt):
        """
        Solves the discretized heat equation implicitely
        with Euler Backward Scheme.
        Input: prev_temp at time it-1
        Returns: new_temp at time it np.array with shape=(nx)
        """
        rcoef = dt / self.rho / self.cp
        cond = self.cond
        dx = self.dx
        prev_temp = np.copy(self.temp)
        matrice = np.zeros((self.nx, self.nx))

        for ix in range(1, self.nx - 1):
            # dkn/dx
            dkn = (cond[ix + 1] - cond[ix]) / (dx[ix]) + (cond[ix] - cond[ix - 1]) / (
                2 * dx[ix - 1]
            )
            # an
            matrice[ix, ix - 1] = (
                -rcoef[ix]
                / dx[ix - 1]
                * (-dkn / 2 + 2 * cond[ix] / (dx[ix] + dx[ix - 1]))
            )
            # bn
            matrice[ix, ix] = 1 - rcoef[ix] / (dx[ix] * dx[ix - 1]) * (
                dkn / 2 * (dx[ix] - dx[ix - 1]) - 2 * cond[ix]
            )
            # cn
            matrice[ix, ix + 1] = (
                -rcoef[ix] / dx[ix] * (dkn / 2 + 2 * cond[ix] / (dx[ix] + dx[ix - 1]))
            )

        # Set BC
        source = prev_temp + rcoef * self.qheat
        matrice, source = self.set_flux_BC(matrice, source, dt)
        # Inverse matrix and compute next iteration
        new_temp = np.linalg.solve(matrice, source)
        return new_temp

    def set_flux_BC(self, matrice, source, dt):
        """
        Set boundary conditions for implicit Euler Scheme
        Imposed flux or imposed temperature possible.
        """
        rcoef = dt / self.rho / self.cp
        cond = self.cond

        # Set Boundary conditions
        bc_top = self.solar_flux / cond[0]
        bc_top += self.eps * cst.SIGMA / self.cond[0] * self.temp[0] ** 4
        self.bc_top = bc_top
        bc_bottom = 0

        source[0] = (
            self.temp[0]
            + rcoef[0] * (cond[1] - 3 * cond[0]) / self.dx[0] * bc_top
            + rcoef[0] * self.qheat[0]
        )
        source[-1] = (
            self.temp[-1]
            + rcoef[-1] * (3 * cond[-1] - cond[-2]) * bc_bottom / self.dx[-1]
            + rcoef[-1] * self.qheat[-1]
        )

        matrice[0, 0] = 1 + 2 * rcoef[0] * cond[0] / self.dx[0] ** 2  # b1
        matrice[0, 1] = -2 * rcoef[0] * cond[0] / self.dx[0] ** 2  # c1
        matrice[-1, -2] = -2 * rcoef[-1] * cond[-1] / self.dx[-1] ** 2  # an
        matrice[-1, -1] = 1 + 2 * rcoef[-1] * cond[-1] / self.dx[-1] ** 2  # bn
        return matrice, source
