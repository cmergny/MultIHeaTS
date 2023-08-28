"""
Run simu when the solar flux has an eclipse,
to see if the stability of the fully implicit scheme
can be advantegeous in these kind of situation in comparison to 
the a semi implicit solver.
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

from multiheats.solvers import ImplicitSolver
from conductionCrankQ import conductionQ, flux_noatm

# from multiheats.visualise import animate_function
from visu import animate_function, animate_2functions


class FakeProfile:
    """Fake profile to reproduce converence figure."""

    def __init__(self):
        """docstring for __init__"""
        self.nx = 100
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
step = 10000
dt = period / step
nt = 5 * step

# Saving dics
dic_times = {}
dic_temps = {}

times = np.zeros(nt)
cn_temps = np.zeros((prof.nx, nt))
im_temps = np.zeros((prof.nx, nt))

prev_temp = np.full(prof.nx, 90.0)  # Schorghofer
prof.temp = np.full(prof.nx, 90.0)  # Me
fluxes = np.zeros(nt)
prev_slr_flux = (1 - prof.alb) * flux_noatm(Rau, decl, lat, HA=0)

for it in tqdm(range(0, nt)):
    # Time
    times[it] = (it + 1) * dt
    HA = 2 * np.pi * ((times[it] / period) % 1.0)
    slr_flux = (1 - prof.alb) * flux_noatm(Rau, decl, lat, HA)
    # Eclipse
    if HA > 0.90 * 2 * np.pi and HA < 2 * np.pi:
        slr_flux = 0

    # Implicit
    solver = ImplicitSolver(prof)
    solver.upBCmethod = "None"
    prof.temp = solver.implicit_scheme(dt, -slr_flux)
    im_temps[:, it] = prof.temp

    cn_temps[:, it] = conductionQ(
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
    fluxes[it] = slr_flux

# Save Everything
# dic_times[step] = times
# dic_temps[step] = temps

# np.save("cn_times", dic_times)
# np.save("cn_temps", dic_temps)

# np.save("im_times", dic_times)
# np.save("im_temps", dic_temps)


eclipse = [int(0.90 * step), step]
start = 0
end = -1
fig, ax = plt.subplots()
ax.plot(times[start:end], fluxes[start:end])
plt.show()

loc = eclipse[0] + 3 * step
idx = loc + 9
fig, ax = plt.subplots()
ax.plot(prof.spaces, cn_temps[:, idx], label="CN")
ax.plot(prof.spaces, im_temps[:, idx], label="Implicit")
plt.legend()
plt.show()

diff_temps = cn_temps - im_temps
fig, ax = plt.subplots()
ax.plot(times, fluxes / fluxes.max(), "--", color="orange", label="Flux")
ax.plot(times, diff_temps[0, :], label=r"$\Delta T$")
ax.set_xlabel("Times")
ax.set_ylabel("Temperature Difference (K)")
plt.legend()
# plt.savefig("eclipse.png")
plt.show()

anim = animate_function(prof.spaces, diff_temps.T, step=5)
plt.show()

diff_temps = cn_temps - im_temps
anim = animate_2functions(
    prof.spaces,
    cn_temps.T,
    im_temps.T,
    "CN",
    "IM",
    step=5,
    # save=True,
    frames=400,
)
plt.show()
