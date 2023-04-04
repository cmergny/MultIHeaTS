""" 
Manually create a surface profile to suit your needs.
(i.e.) spaces, times, density, thermic capacity, conductivity, heat source.
"""

### IMPORTS
import numpy as np


class Profile:
    """
    Surface profile used by the termal model.
    You may change properties directly in the code.

    """

    def __init__(self) -> None:
        self.nx = 100
        self.lat = 0
        self.long = 0
        self.eps = 0.94  # Emissivity
        x0 = 0  # Surface depth (m)
        xf = 5  # Total depth (m)

        self.qheat = np.full(self.nx, 0)
        power = 3
        spaces = np.linspace(x0, xf ** (1 / power), self.nx)
        self.spaces = spaces ** (power)
        self.spaces = np.linspace(0, xf, self.nx)

    def monolayer_prof(self):
        """
        Generate an monolayered surface profile.
        PARAMS:
            cond - Conductivity (W.m-1.K-2)
            rho - Density (kg,m-3)
            cp - Heat capacity (J.kg-1.K-1)
        """
        cond = 0.03
        rho = 917.0
        cp = 839.0

        self.cond = np.full(self.nx, cond)
        self.rho = np.full(self.nx, rho)
        self.cp = np.full(self.nx, cp)
        self.interf = 0

    def bilayer_prof(self):
        """
        Generate an bilayered surface profile.
        PARAMS:
            interface - in meters, controls the thickness of the top layer. (must be > xf)
            transition - thickness of the smooth transition between the top and bottom layer.
        """
        # TOP
        cond_top = 0.01
        rho_top = 917.0
        cp_top = 839.0
        # BOTTOM
        cond_bot = cond_top / 2
        rho_bot = rho_top / 2
        cp_bot = cp_top / 2

        thermal_skin = self.thermal_skin(cond_top, rho_top, cp_top)
        interface = 2 * thermal_skin  # (m)
        transition = 2e-3  # (m)

        self.cond = bilayer(self.spaces, cond_top, cond_bot, interface, transition)
        self.rho = bilayer(self.spaces, rho_top, rho_bot, interface, transition)
        self.cp = bilayer(self.spaces, cp_top, cp_bot, interface, transition)
        self.interf = interface

    def thermal_skin(self, cond, rho, cp):
        """delta = (2*alpha/omega)**1/2"""
        period_day = 79.33 * 86400  # Japet Day
        omega = 2 * np.pi / period_day
        alpha = cond / rho / cp
        self.skin = np.sqrt(2 * alpha / omega)
        return self.skin

    def stability_number(self, dt):
        """Compute r = alpha dt/dx^2"""
        dx = np.diff(self.spaces).min()
        alpha = np.min(self.cond / self.rho / self.cp)
        return alpha * dt / dx**2


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
