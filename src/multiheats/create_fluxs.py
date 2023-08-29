import numpy as np

import multiheats.constants as cst


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
    solar_flux[np.where(solar_flux > 0)] = 0
    return solar_flux


def get_eq_temp(solar_fluxs, alb, eps):
    """Rough Estimate of Initial Temperature"""
    temp_eq = (solar_fluxs.mean() / eps / cst.SIGMA) ** (1 / 4)
    return temp_eq
