"""
Implicit solver of the heat equation.
"""

### IMPORTS
import numpy as np

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
        # new_temp = np.linalg.solve(matrice, source)
        matric_inv = np.linalg.inv(matrice)
        new_temp = np.dot(matric_inv, source)

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


class OldImplicitSolver(ImplicitSolver):
    """Old class"""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.copy_temp = np.copy(self.temp)

    def implicit_scheme(self, nt, dt):
        """
        Solves the discretized heat equation implicitely with Euler Backward Scheme.
        Returns: u[:, :] np.array with u.shape=(nx, nt)
        """
        rcoef = dt / self.rho / self.cp
        self.rcoef = rcoef
        A = np.zeros((self.nx, self.nx, nt))
        self.temp = np.zeros((self.nx, nt))
        self.temp[:, 0] = self.copy_temp

        for it in range(1, nt):
            for ix in range(1, self.nx - 1):
                A[ix, ix - 1] = (
                    rcoef[ix]
                    / self.dx[ix] ** 2
                    / 4
                    * (self.cond[ix + 1] - self.cond[ix - 1] - 4 * self.cond[ix])
                )  # an
                A[ix, ix] = 1 + 2 * rcoef[ix] / self.dx[ix] ** 2 * self.cond[ix]  # bn
                A[ix, ix + 1] = (
                    rcoef[ix]
                    / 4
                    * self.dx[ix] ** 2
                    * (self.cond[ix - 1] - self.cond[ix + 1] - 4 * self.cond[ix])
                )  # cn
            # Add source term
            source = np.copy(self.temp[:, it - 1]) + rcoef[ix] * 0
            # Set BC
            A, source = self.set_BC(A, source, it)
            # Inverse matrix
            Ainv = np.linalg.inv(A[:, :, it])
            self.temp[:, it] = np.dot(Ainv, source)
        return self.temp

    def set_BC(self, A, source, it, mode="Neumman"):
        """Set boundary conditions for implicit Euler Scheme"""
        if mode == "Neumman":
            source[0] = self.temp[0, it - 1]
            source[-1] = self.temp[-1, it - 1]

            # Source = T-1 si BC=0 et Q=0
            # Matrix
            A[0, 0, :] = 1 + 2 * self.rcoef[0] * self.cond[0] / self.dx[0] ** 2  # b1
            A[0, 1, :] = -2 * self.rcoef[0] * self.cond[0] / self.dx[0] ** 2  # c1
            A[self.nx - 1, self.nx - 2, :] = (
                -2 * self.rcoef[-1] * self.cond[-1] / self.dx[-1] ** 2
            )  # an
            A[self.nx - 1, self.nx - 1, :] = (
                1 + 2 * self.rcoef[-1] * self.cond[-1] / self.dx[-1] ** 2
            )  # bn
        # Change mode
        elif mode == "Dirichlet":
            source[0] = self.BC[0, it]
            source[-1] = self.BC[1, it]
            # Matrix
            A[0, 0, :] = 1
            A[0, 1, :] = 0
            A[self.nx - 1, self.nx - 2, :] = 0
            A[self.nx - 1, self.nx - 1, :] = 1
        return (A, source)
