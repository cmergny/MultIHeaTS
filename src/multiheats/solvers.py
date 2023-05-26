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
        self.need_update = True  # Update the matrix
        self.upBCmethod = "Old"  # Correction Upper BC

    def implicit_scheme(self, dt, solar_flux):
        """
        Solves the discretized heat equation implicitely
        with Euler Backward Scheme.
        Input: prev_temp at time it-1
        Returns: new_temp at time it np.array with shape=(nx)
        """

        rcoef = dt / self.rho / self.cp

        if self.need_update:
            self.mat = self.update_matrix(self.nx, self.cond, self.dx, rcoef)
            self.inv_mat = np.linalg.inv(self.mat)

        # Set BC
        source = self.temp + rcoef * self.qheat

        s0, sN = self.set_flux_BC(rcoef, solar_flux)
        source[0] = s0
        source[-1] = sN

        new_temp = np.matmul(self.inv_mat, source)

        return new_temp

    def set_flux_BC(self, rcoef, solar_flux):
        """
        Set boundary conditions for implicit Euler Scheme
        Imposed flux or imposed temperature possible.
        """

        self.bc_bottom = 0
        s0_corr, b0_corr = self.apply_upBCmethod(
            self.upBCmethod, rcoef
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

        return s0, sN

    def apply_upBCmethod(self, method, rcoef):
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

        # BC values
        coef_top = rcoef[0] * self.cond[0] / self.dx[0] ** 2
        coef_bot = rcoef[-1] * self.cond[-1] / self.dx[-1] ** 2
        b0 = 1 + 2 * coef_top  # + b0_corr  # b0
        c0 = -2 * coef_top  # c0
        aN = -2 * coef_bot  # aN
        bN = 1 + 2 * coef_bot  # bN
        # Apply to matric
        matrix[1, 0] = b0
        matrix[1, -1] = bN
        matrix[0, 1] = c0
        matrix[2, -2] = aN

        # Real matrice
        mat = (
            np.diag(matrix[2, :-1], k=-1)
            + np.diag(matrix[1])
            + np.diag(matrix[0, 1:], k=1)
        )
        return mat
