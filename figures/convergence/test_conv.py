import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

from conductionCrankQ import conductionQ, flux_noatm


# ext_times = np.linspace(0, tf + dt, nt + 1)
# temps = np.zeros((prof.nx, nt))
# temps[:, 0] = prof.temp


class FakeProfile:
    """Fake profile to reproduce converence figure."""

    def __init__(self):
        """docstring for __init__"""
        self.nx = 31
        self.spaces = np.load("depth.npy")
        self.inertia = np.full(self.nx, 120)
        self.eps = 1


# Params
prof = FakeProfile()
period = 88775.244 * 670
prof.rho_cp = prof.inertia * np.sqrt(period / np.pi)
# Flux
lat = 5 * np.pi / 180
Rau = 1.52
decl = 0

# Steps
coef = np.arange(7, 8)
steps = 2 ** (12 - coef)
dts = period / steps
# Time
nts = 512 * steps + 1  # Add one step bc we count init temp

for isp, step in enumerate(steps):
    nt = nts[isp]
    dt = dts[isp]
    times = np.zeros(nt)
    temps = np.zeros((prof.nx, nt))
    temps[:, 0] = 210

    for it in tqdm(range(1, nt)):
        # time
        times[it] = it * dt
        # Use prev HA ?
        HA = 2 * np.pi * ((times[it] / period) % 1)
        slr_flux = flux_noatm(Rau, decl, lat, HA)
        next_HA = 2 * np.pi * (((times[it] + dt) / period) % 1)
        next_slr_flux = flux_noatm(Rau, decl, lat, next_HA)

        temps[:, it] = conductionQ(
            prof.nx - 1,
            prof.spaces,
            dt,
            slr_flux,
            next_slr_flux,
            temps[:, it - 1],
            prof.inertia,
            prof.rho_cp,
            prof.eps,
            0,
            0,
        )

    print(f"time (days) = {times[-1]/86400}")


idx = -1
fig, ax = plt.subplots()
ax.plot(prof.spaces, temps[:, idx])
plt.show()
