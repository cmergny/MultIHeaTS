import numpy as np
import multiheats.constants as cst


def guess_ini_temp(period, dt, distance, diurn_per, alb, eps):
    """
    Compute Equilibrium Temperature after a given period
    of time and a timestep. Based on the equilibrium of
    energy bw solar flux and black body emission at the surface.
    RETURNS:
        temp_eq: float, equilibrium temperature
    """
    times = np.arange(0, period, dt)  # All iters during the period
    solar_flux = (
        (1 - alb)
        * cst.SOLAR_CST
        / (distance / cst.UA) ** 2
        * np.cos(2 * np.pi * times / diurn_per + np.pi)
    )
    solar_flux[np.where(solar_flux > 0)] = 0
    temp_eq = (-solar_flux.mean() / eps / cst.SIGMA) ** (1 / 4)
    return temp_eq


def fake_slr_flux(alb, times, distance, period):
    """
    Use a truncated cosinus function to modelise the solar flux.
    alb - albedo
    times - time array (in s)
    distance - distance from sun (m)
    period - diurnal period (s)
    Note: flux is negative because counted from the surface to the exterior.
    """
    solar_flux = (
        (1 - alb)
        * cst.SOLAR_CST
        / (distance / cst.UA) ** 2
        * np.cos(2 * np.pi * times / period + np.pi)
    )
    if solar_flux > 0:
        solar_flux = 0
    return -solar_flux  # negative sign
