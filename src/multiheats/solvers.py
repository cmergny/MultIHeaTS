"""
Implicit solver of the heat equation.
"""

### IMPORTS
import numpy as np
import scipy

import multiheats.constants as cst

### CLASS


class ImplicitSolver:
    """
    Solver class which uses an implicit
    method to solve the heat equation.
    """

    def __init__(self, prof) -> None:
        self.temp = prof.temp
        self.cond = prof.cond
        self.rho = prof.rho
        self.cp = prof.cp
        self.qheat = prof.qheat
        self.eps = prof.eps
        self.nx = prof.spaces.shape[0]
        self.dx = np.diff(prof.spaces)
        self.need_update = True

    def implicit_scheme(self, dt, solar_flux):
        """
        Solves the discretized heat equation implicitely
        with Euler Backward Scheme.
        Input: prev_temp at time it-1
        Returns: new_temp at time it np.array with shape=(nx)
        """

        rcoef = dt / self.rho / self.cp

        if self.need_update:
            self.matrix = self.update_matrix(self.nx, self.cond, self.dx, rcoef)

        # Set BC
        source = self.temp + rcoef * self.qheat

        s1, sN, b1, c1, aN, bN = self.set_flux_BC(rcoef, solar_flux)
        source[0] = s1
        source[-1] = sN
        self.matrix[1, 0] = b1
        self.matrix[1, -1] = bN
        self.matrix[0, 1] = c1
        self.matrix[2, -2] = aN

        new_temp = scipy.linalg.solve_banded(
            (1, 1),
            self.matrix,
            source,
            check_finite=False,
            overwrite_b=True,
        )

        return new_temp

    def set_flux_BC(self, rcoef, solar_flux):
        """
        Set boundary conditions for implicit Euler Scheme
        Imposed flux or imposed temperature possible.
        """

        # Set Boundary conditions
        self.bc_top = (
            solar_flux + self.eps * cst.SIGMA * self.temp[0] ** 4
        ) / self.cond[0]

        self.bc_bottom = 0

        s1 = self.temp[0] + rcoef[0] * (
            (self.cond[1] - 3 * self.cond[0]) / self.dx[0] * self.bc_top + self.qheat[0]
        )

        sN = self.temp[-1] + rcoef[-1] * (
            (3 * self.cond[-1] - self.cond[-2]) * self.bc_bottom / self.dx[-1]
            + self.qheat[-1]
        )

        coef_top = rcoef[0] * self.cond[0] / self.dx[0] ** 2
        coef_bot = rcoef[-1] * self.cond[-1] / self.dx[-1] ** 2

        b1 = 1 + 2 * coef_top  # b1
        c1 = -2 * coef_top  # c1
        aN = -2 * coef_bot  # aN
        bN = 1 + 2 * coef_bot  # bN
        return s1, sN, b1, c1, aN, bN

    def update_matrix(self, nx, cond, dx, rcoef):
        """
        Update the coefficients of the tridiagonal matrix.
        Necessary if some of the thermal properties, the grid or the timestep are changed between iterations.
        """
        matrix = np.zeros((3, nx))
        dkn = (cond[2:] - cond[1:-1]) / (2 * dx[1:]) + (cond[1:-1] - cond[:-2]) / (
            2 * dx[:-1]
        )
        # an
        matrix[2, :-2] = (
            -rcoef[1:-1] / dx[:-1] * (-dkn / 2 + 2 * cond[1:-1] / (dx[1:] + dx[:-1]))
        )
        # bn
        matrix[1, 1:-1] = 1 - rcoef[1:-1] / (dx[1:] * dx[:-1]) * (
            dkn / 2 * (dx[1:] - dx[:-1]) - 2 * cond[1:-1]
        )
        # cn
        matrix[0, 2:] = (
            -rcoef[1:-1] / dx[1:] * (dkn / 2 + 2 * cond[1:-1] / (dx[1:] + dx[:-1]))
        )
        return matrix
