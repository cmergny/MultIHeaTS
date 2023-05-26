"""
Implicit solver of the heat equation.
"""

### IMPORTS
import cupy as np

# import cupyx.scipy as scipy
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
        self.need_update = True  # Update the matrix
        self.upBCmethod = "Leyrat"  # Correction Upper BC

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

        s0, sN, b0, c0, aN, bN = self.set_flux_BC(rcoef, solar_flux)
        source[0] = s0
        source[-1] = sN
        self.matrix[1, 0] = b0
        self.matrix[1, -1] = bN
        self.matrix[0, 1] = c0
        self.matrix[2, -2] = aN

        # new_temp = scipy.linalg.solve_banded(
        #     (1, 1),
        #     self.matrix,
        #     source,
        #     check_finite=False,
        #     overwrite_b=True,
        # )
        new_temp = self.temp

        return new_temp

    def set_flux_BC(self, rcoef, solar_flux):
        """
        Set boundary conditions for implicit Euler Scheme
        Imposed flux or imposed temperature possible.
        """

        self.bc_bottom = 0
        s0_corr, b0_corr = self.apply_upBCmethod(
            self.upBCmethod, rcoef, solar_flux
        )  # Correction terms

        s0 = self.temp[0] + rcoef[0] * (
            (self.cond[1] - 3 * self.cond[0])
            / self.dx[0]
            * (solar_flux / self.cond[0] + s0_corr)
            + self.qheat[0]
        )

        sN = self.temp[-1] + rcoef[-1] * (
            (3 * self.cond[-1] - self.cond[-2]) * self.bc_bottom / self.dx[-1]
            + self.qheat[-1]
        )

        coef_top = rcoef[0] * self.cond[0] / self.dx[0] ** 2
        coef_bot = rcoef[-1] * self.cond[-1] / self.dx[-1] ** 2

        b0 = 1 + 2 * coef_top + b0_corr  # b0
        c0 = -2 * coef_top  # c0
        aN = -2 * coef_bot  # aN
        bN = 1 + 2 * coef_bot  # bN
        return s0, sN, b0, c0, aN, bN

    def apply_upBCmethod(self, method, rcoef, solar_flux):
        """
        Add corrected terms for the upper boundary conditions.
        Following C. Leyrat's solver (Probably same as Schorghofer).
        If not uses our standard BC (Mergny et al 2023), may cause instabilities.
        """
        if method == "Leyrat":
            s0_corr = (
                -3 * cst.SIGMA * self.eps * self.temp[0] ** 4 / self.cond[0]
            )  # Correction from Leyrat BC
            b0_corr = (
                -4
                * rcoef[0]
                * (self.cond[1] - 3 * self.cond[0])
                * cst.SIGMA
                * self.eps
                * self.temp[0] ** 3
                / (self.dx[0] * self.cond[0])
            )  # Correction with Leyrat BC
        else:
            s0_corr = self.eps * cst.SIGMA * self.temp[0] ** 4 / self.cond[0]
            b0_corr = 0
        return s0_corr, b0_corr

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
