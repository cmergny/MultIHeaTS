import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

from multiheats.solvers import ImplicitSolver
from conductionCrankQ import conductionQ, flux_noatm


class FakeProfile:
    """Fake profile to reproduce converence figure."""

    def __init__(self):
        """docstring for __init__"""
        self.nx = 101
        self.length = 2
        # Remember that z[2]=3*z[1]
        self.spaces = np.linspace(0, self.length, self.nx)
        self.spaces[1] = self.spaces[2] / 3
        self.inertia = np.full(self.nx, 200)
        self.eps = 1
        self.alb = 0.015


# Params
prof = FakeProfile()
period = 79.3 * 86400
prof.rho = np.full(prof.nx, 800.0)
prof.cp = np.full(prof.nx, 600.0)
prof.cond = prof.inertia**2 / prof.rho / prof.cp
prof.qheat = np.zeros(prof.nx)
prof.rho_cp = prof.rho * prof.cp

# Flux
lat = 5 * np.pi / 180
Rau = 9.51
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

    # prev_temp = np.full(prof.nx, 90.0)  # Schorghofer
    prof.temp = np.full(prof.nx, 90.0)  # Me
    prev_slr_flux = (1 - prof.alb) * flux_noatm(Rau, decl, lat, HA=0)

    for it in tqdm(range(0, nt)):
        # Time
        times[it] = (it + 1) * dt
        HA = 2 * np.pi * ((times[it] / period) % 1.0)
        slr_flux = (1 - prof.alb) * flux_noatm(Rau, decl, lat, HA)

        solver = ImplicitSolver(prof)
        prof.temp = solver.implicit_scheme(dt, -slr_flux)
        temps[:, it] = prof.temp

        # temps[:, it] = conductionQ(
        #     prof.nx - 1,
        #     prof.spaces,
        #     dt,
        #     prev_slr_flux,
        #     slr_flux,
        #     prev_temp,
        #     prof.inertia,
        #     prof.rho_cp,
        #     prof.eps,
        #     0,
        #     0,
        # )
        # prev_slr_flux = slr_flux

    # Save Everything
    dic_times[step] = times
    dic_temps[step] = temps

# np.save("cn_times", dic_times)
# np.save("cn_temps", dic_temps)

# np.save("im_times", dic_times)
# np.save("im_temps", dic_temps)

# idx = -1
# fig, ax = plt.subplots()
# ax.plot(prof.spaces, temps[:, idx])
# plt.show()
