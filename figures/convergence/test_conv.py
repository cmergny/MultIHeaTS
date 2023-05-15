import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

from conductionCrankQ import conductionQ, flux_noatm


class FakeProfile:
    """Fake profile to reproduce converence figure."""

    def __init__(self):
        """docstring for __init__"""
        self.nx = 31
        self.spaces = np.load("depth.npy")
        self.inertia = np.full(self.nx, 120)
        self.eps = 1
        self.alb = 0.2


# Params
prof = FakeProfile()
period = 88775.244 * 670
prof.rho_cp = prof.inertia * np.sqrt(period / np.pi)
# Flux
lat = 5 * np.pi / 180
Rau = 1.52
decl = 0

# Steps
coef = np.arange(0, 8)
steps = 2 ** (12 - coef)
dts = period / steps
# Time
nts = 512 * steps  # Add one step bc we count init temp

# Saving dics
dic_times = {}
dic_temps = {}

for isp, step in enumerate(steps):
    nt = nts[isp]
    dt = dts[isp]
    times = np.zeros(nt)
    temps = np.zeros((prof.nx, nt))

    prev_temp = np.full(prof.nx, 210.0)
    prev_slr_flux = (1 - prof.alb) * flux_noatm(Rau, decl, lat, HA=0)

    for it in tqdm(range(0, nt)):
        # Time
        times[it] = (it + 1) * dt
        HA = 2 * np.pi * ((times[it] / period) % 1.0)
        slr_flux = (1 - prof.alb) * flux_noatm(Rau, decl, lat, HA)

        temps[:, it] = conductionQ(
            prof.nx - 1,
            prof.spaces,
            dt,
            prev_slr_flux,
            slr_flux,
            prev_temp,
            prof.inertia,
            prof.rho_cp,
            prof.eps,
            0,
            0,
        )
        prev_slr_flux = slr_flux

    # Save Everything
    dic_times[step] = times
    dic_temps[step] = temps

np.save("dic_times", dic_times)
np.save("dic_temps", dic_temps)

# idx = -1
# fig, ax = plt.subplots()
# ax.plot(prof.spaces, temps[:, idx])
# plt.show()
