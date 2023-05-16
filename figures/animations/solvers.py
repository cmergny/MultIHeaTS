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


class CrankNicolson:
    """
    Norbert Schorghofer Implementation of CN of 1D heat equation.
    Modified by C. Mergny 2023 to fit into the multiheats pipeline.
    originally written by Samar Khatiwala, 2001
    extended to variable thermal properties
    and irregular grid by Norbert Schorghofer, 2004
    added predictor-corrector 9/2019
    converted to Python 3/2021
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
        self.spaces = prof.spaces

    def CN_scheme(self, dt, solar_flux, next_solar_flux):
        """
        conductionQ:  program to calculate the diffusion of temperature into the ground and thermal emission at the surface with variable thermal properties on irregular grid Crank-Nicolson scheme, flux conservative uses Samar's radiation formula.

        Eqn: rhoc*T_t = (k*T_z)_z
        BC (z=0): Q(t) + kT_z = em*sig*T^4
        BC (z=L): heat flux = Fgeotherm

        nx = number of grid points (not counting the surface)
        dt = time step
        solar_flux,next_solar_flux = net solar insolation at time steps n and n+1 [W/m^2] # the albedo has already been applied.
        T = vertical temperature profile [K]  (in- and output)
        T[0] = surface temperature [K]  (in- and output)
        ti = thermal inertia [J m^-2 K^-1 s^-1/2]  VECTOR
        rhoc = rho*c  VECTOR where rho=density [kg/m^3] and c=specific heat [J K^-1 kg^-1]
        ti and rhoc are not allowed to vary in the layers immediately adjacent to the surface or the bottom
        Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]
        Fsurf = heat flux at surface [W/m^2]  (output)

        Grid: surface is at z=0
              z[0]=0, z[2]=3*z[1], i.e., the top layer has half the width
              T[1] is at z[1]; ...; T[i] is at z[i]
              k[i], rhoc[i], ti[i] are midway between z[i-1] and z[i]

        """
        rhoc = self.rho * self.cp
        nx = self.nx
        bc_bottom = 0

        # set some constants
        dz = self.spaces[1]
        alpha = np.zeros(nx)
        gamma = np.zeros(nx)
        beta = dt / (rhoc[1] + rhoc[2]) / dz**2
        alpha[1] = beta * self.cond[2]
        gamma[1] = beta * self.cond[1]

        for i in range(2, nx - 1):
            buf = dt / (self.spaces[i + 1] - self.spaces[i - 1])
            alpha[i] = (
                2
                * self.cond[i + 1]
                * buf
                / (rhoc[i] + rhoc[i + 1])
                / (self.spaces[i + 1] - self.spaces[i])
            )
            gamma[i] = (
                2
                * self.cond[i]
                * buf
                / (rhoc[i] + rhoc[i + 1])
                / (self.spaces[i] - self.spaces[i - 1])
            )
        buf = dt / (self.spaces[-1] - self.spaces[-2]) ** 2
        gamma[-1] = self.cond[-1] * buf / (2 * rhoc[-1])  # assumes rhoc[nx+1]=rhoc[nx]
        alpha[-1] = 0.0  # ensures b[nx-1] = 1 + gamma[nx-1]

        self.cond1 = self.cond[0] / dz

        # nx-1 because T[0] traité à part.
        matrice = np.zeros((3, nx - 1))
        matrice[0, 1:] = -alpha[1:-1]  # coefficient 'c'
        matrice[1, :] = 1.0 + alpha[1:] + gamma[1:]  # coefficient 'b'
        matrice[2, :-1] = -gamma[2:]  # coefficient 'a'

        Tr = self.temp[0]  # 'reference' temperature
        itr = 0
        new_temp = np.copy(self.temp)

        while itr < 10:

            if itr > 0:
                Tr = np.sqrt(
                    Tr * new_temp[0]
                )  # lineariself.spacese around an intermediate temperature
                new_temp[1:] = self.temp[1:]
                print(f"Reiterating {itr}")

            # Emission
            arad = -3.0 * self.eps * cst.SIGMA * Tr**4
            brad = 2.0 * self.eps * cst.SIGMA * Tr**3
            ann = (solar_flux - arad) / (self.cond1 + brad)
            annp1 = (next_solar_flux - arad) / (self.cond1 + brad)
            bn = (self.cond1 - brad) / (self.cond1 + brad)
            b1 = 1.0 + alpha[1] + gamma[1] - gamma[1] * bn  # b[1]

            # Set RHS
            source = np.zeros(nx)
            source[1] = (
                gamma[1] * (annp1 + ann)
                + (1.0 - alpha[1] - gamma[1] + gamma[1] * bn) * new_temp[1]
                + alpha[1] * new_temp[2]
            )
            for i in range(2, nx - 1):  # 2...nx-1
                source[i] = (
                    gamma[i] * new_temp[i - 1]
                    + (1.0 - alpha[i] - gamma[i]) * new_temp[i]
                    + alpha[i] * new_temp[i + 1]
                )
            source[-1] = (
                gamma[-1] * new_temp[-2]
                + (1.0 - gamma[-1]) * new_temp[-1]
                + dt / rhoc[-1] * bc_bottom / (self.spaces[-1] - self.spaces[-2])
            )  # assumes rhoc[nx+1]=rhoc[nx]

            matrice[1, 0] = b1  # coefficient b[1]

            # Solve for T at n+1
            new_temp[1:] = scipy.linalg.solve_banded((1, 1), matrice, source[1:])
            new_temp[0] = 0.5 * (annp1 + bn * new_temp[1] + new_temp[1])  # (T0+T1)/2

            # iterative predictor-corrector
            if new_temp[0] < 1.2 * Tr and new_temp[0] > 0.8 * Tr:
                break
            itr += 1

        if itr > 9:
            print("Did not find convergence")

        return new_temp
