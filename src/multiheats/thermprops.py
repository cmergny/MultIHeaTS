# Compute thermal properties of porous ice
import numpy as np


def get_cond_bulk(temp):
    """
    Conductivity of bulk ice with temperature, Ferrari et al 2016
    Temperature must be negative.
    """
    return 567 / temp


def get_cond_porous(rg, rb, phi, cond_bulk):
    """
    Gundlach2011b, from Chan Tien 1973, only valid when rb<<rg?
    Using Eq 12 where the division by rg comes from Eps*rg in Eq.18
    Hertz expression for rb does not need to be computed again
    """
    # Epsilon Eq. 18 Gund
    f1 = 5.18e-2  # w/o uncertainties
    f2 = 5.26  # w/o uncertainties
    epsilon_r = f1 * np.exp(f2 * (1 - phi))  # Epsilon*rg = adim
    return cond_bulk * epsilon_r * rb / rg


def old_get_cond(rg, rb, phi, cond_bulk):
    """
    Conductivity of porous ice (Mergny2024, LunaIcy).
    rb must be < to rg.
    """
    return cond_bulk * (1 - phi) * rb / rg


def thermal_skin(cond, rho, cp, period):
    """delta = (2*alpha/omega)**1/2"""
    omega = 2 * np.pi / period
    alpha = cond / rho / cp
    return np.sqrt(2 * alpha / omega)


def ice_density_depth(x, phi0, gravity):
    """
    Compute the density profile as a function of depth (Mergny et al. 2024)
    lambda (lbd) parameter obtain through fit of function (Wilson 1994, eq4)
    used on polar ice data (Horhold 1994, fig 8)
    """
    lbd = 4.918407714021047e-06
    rho_bulk = 917  # density pure ice SI (valid for T=0 C only!!)
    rho = rho_bulk / (1 + (phi0 / (1 - phi0)) * np.exp(-lbd * rho_bulk * gravity * x))
    return rho


if __name__ == "__main__":
    phis = 0.2  # Porosity
    temp = 100  # temperature (K)
    rg = 100e-6  # grain radius (m)
    rb = 40e-6  # neck radius (m)

    cond_bulk = get_cond_bulk(temp)
    cond_me = get_cond_porous(rg, rb, phis, cond_bulk)
    print(cond_bulk, cond_me)
