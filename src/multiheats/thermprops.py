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
    Conductivity of porous ice (Mergny2024, LunaIcy).
    rb must be < to rg.
    """
    cond = cond_bulk * (1 - phi) * rb / rg

    return cond


def ice_density_depth(x, phi0, gravity):
    """
    Compute the density profile as a function of depth (Mergny et al. 2024)
    lambda (lbd) parameter obtain through fit of function (Wilson 1994, eq4)
    used on polar ice data (Horhold 1994, fig 8)
    """
    lbd = 4.918407714021047e-06
    rho_bulk = 917  # density pure ice SI (valid for T=0 C only!!)
    rho = rho_bulk / (1 + (phi0 / (1 - phi0)) *
                      np.exp(-lbd * rho_bulk * gravity * x))
    return rho


if __name__ == "__main__":
    phis = 0.2  # Porosity
    temp = 100  # temperature (K)
    rg = 100e-6  # grain radius (m)
    rb = 40e-6  # neck radius (m)

    cond_bulk = get_cond_bulk(temp)
    cond_me = get_cond_porous(rg, rb, phis, cond_bulk)
    print(cond_bulk, cond_me)
