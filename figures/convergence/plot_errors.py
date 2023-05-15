import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

dic_times = np.load("dic_times.npy", allow_pickle=True).item()
dic_temps = np.load("dic_temps.npy", allow_pickle=True).item()

# Get params
steps = np.array(list(dic_times.keys()))
ns = steps.size
dts = np.zeros(ns)
for isp, val in enumerate(dic_times.values()):
    dts[isp] = val[0]

# Get error
dic_errs = {}
for isp, step in enumerate(steps[1:]):
    prev_step = steps[isp]  # isp starts at 0
    temps = dic_temps[step]
    prev_temps = dic_temps[prev_step]
    times = dic_times[step]
    prev_times = dic_times[prev_step]

    # Check that times are coherent
    err_time = abs(prev_times[1::2] - times)
    print("Time diff ", err_time, err_time.max(), err_time.min())

    # Store Temperature Diff
    dic_errs[step] = prev_temps[:, 1::2] - temps

err_temps = np.zeros(ns - 1)
for i, val in enumerate(dic_errs.values()):
    err_temps[i] = np.abs(val)[:, -1].mean()


period = 88775.244 * 670
dt_frac = dts / period

fig, ax = plt.subplots()
ax.plot(dt_frac[1:], err_temps, ".", label="Mean Error")
ax.plot(dt_frac, 5e2 * dt_frac**2, label="dt^2")
ax.plot(dt_frac, 5e2 * dt_frac, label="dt")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("dt")
ax.set_ylabel("Error in K")
# ax.set_ylim([1e-4, 1e2])
ax.set_xlim([1e-4, 1e-1])
plt.legend()
plt.show()
