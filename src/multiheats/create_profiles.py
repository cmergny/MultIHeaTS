"""
Manually create a surface profile to suit your needs.
(i.e.) spaces, times, density, thermic capacity, conductivity, heat source.
"""

# IMPORTS
import numpy as np


class Profile:
    """
    Surface profile used by the termal model.
    You may change properties directly in the code.

    """

    def __init__(self, nx, eps, xmin, xmax, power=1, lat=None, long=None):
        """
        Initialize profile
        lat/long only necessary when importing data
        """
        self.nx = nx  # Nbr of grid points
        self.eps = eps  # Emissivity
        # Optional
        self.qheat = np.full(self.nx, 0)  # Heat source terms
        self.long = long  # Longitude
        # Even or Uneven depths array
        self.spaces = np.linspace(
            xmin, xmax ** (1 / power), self.nx) ** (power)

    def monolayer_prof(self, cond, rho, cp):
        """
        Generate an monolayered surface profile.
        PARAMS:
            cond - Conductivity (W.m-1.K-2)
            rho - Density (kg,m-3)
            cp - Heat capacity (J.kg-1.K-1)
        """
        self.cond = np.full(self.nx, cond)
        self.rho = np.full(self.nx, rho)
        self.cp = np.full(self.nx, cp)
        self.interf = 0

    def bilayer_prof(
        self,
        cond_top,
        rho_top,
        cp_top,
        cond_bot,
        rho_bot,
        cp_bot,
        interface,
        transition=2e-3,
    ):
        """
        Generate an bilayered surface profile.
        PARAMS:
            interface - in meters, controls the thickness of the top layer. (must be > xf)
            transition - thickness of the smooth transition between the top and bottom layer.
        """
        self.cond = bilayer(self.spaces, cond_top,
                            cond_bot, interface, transition)
        self.rho = bilayer(self.spaces, rho_top, rho_bot,
                           interface, transition)
        self.cp = bilayer(self.spaces, cp_top, cp_bot, interface, transition)
        self.interf = interface

    def thermal_skin(self, period):
        """delta = (2*alpha/omega)**1/2"""
        omega = 2 * np.pi / period
        alpha = self.cond / self.rho / self.cp
        self.skin = np.sqrt(2 * alpha / omega)
        return self.skin

    def stability_number(self, dt):
        """Compute r = alpha dt/dx^2"""
        dx = np.diff(self.spaces).min()
        alpha = np.min(self.cond / self.rho / self.cp)
        return alpha * dt / dx**2

    def compute_inertia(self):
        return np.sqrt(self.cond * self.cp * self.rho)


def bilayer(x, start, end, inter, transition=1):
    """
    Returns a bilayer profile (array) with a smooth
    transition (tanh) between a 'start' and 'end' values
    located depth 'inter' and of thickness 'transition'.
    PARAMS:
        x - array representing the space discretisation.
        transition - in meters, thickness of the transition bw top and bot.
    """
    xprime = (x - inter) / x.max() * 1 / transition
    res = (start + end) / 2 + (end - start) / 2 * np.tanh(xprime)
    return res
