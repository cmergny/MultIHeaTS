import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt


def err_with_previous(dic_temps, dic_times):
    """Compute the difference of temperature."""
    dic_errs = {}
    for isp, step in enumerate(steps[1:]):
        prev_step = steps[isp]  # isp starts at 0
        temps = dic_temps[step]
        prev_temps = dic_temps[prev_step]
        # Store Temperature Diff
        dic_errs[step] = prev_temps[:, 1::2] - temps

        # Check that times are coherent
        times = dic_times[step]
        prev_times = dic_times[prev_step]
        err_time = abs(prev_times[1::2] - times)
        print("Time diff ", err_time, err_time.max(), err_time.min())

    return dic_errs


def err_with_ref(dic_temps, dic_times, ref_temps, ref_times):
    """Compute the difference of temperature with a reference temp."""
    dic_errs = {}
    for step in steps[1:]:
        # Temperature Diff
        temps = dic_temps[step]
        jump = int(4096 / step)  # jump bw the ref
        dic_errs[step] = ref_temps[:, jump - 1 :: jump] - temps

        # Check that times are coherent
        times = dic_times[step]
        print(jump)
        err_time = abs(ref_times[jump - 1 :: jump] - times)
        print("Time diff ", err_time, err_time.max(), err_time.min())

    return dic_errs


def get_err_it(dic_errs, it):
    """Return the error for the iteration it"""
    ns = len(dic_errs.keys()) + 1
    err_temps = np.zeros(ns - 1)
    for i, val in enumerate(dic_errs.values()):
        err_temps[i] = np.abs(val)[:, it].mean()
    return err_temps


def get_err_max(dic_errs):
    """Return the error for the iteration it"""
    ns = len(dic_errs.keys()) + 1
    err_temps = np.zeros(ns - 1)
    for i, val in enumerate(dic_errs.values()):
        # err_temps[i] = np.abs(val)[:, :].max()
        err_temps[i] = np.max(np.abs(val)[:, :], axis=0).mean()
    return err_temps


def get_err_mean(dic_errs):
    """Return the error for the iteration it"""
    ns = len(dic_errs.keys()) + 1
    err_temps = np.zeros(ns - 1)
    for i, val in enumerate(dic_errs.values()):
        err_temps[i] = np.abs(val)[:, :].mean()
    return err_temps


def repro_conv_scho(dic_errs):
    """Reproduces the plot that Schorghofer did by email."""
    err_temps = get_err_it(dic_errs, it=-1)

    fig, ax = plt.subplots()
    ax.plot(dt_frac[1:], err_temps, ".", label="Mean Error")
    ax.plot(dt_frac, 5e2 * dt_frac**2, label="dt^2")
    ax.plot(dt_frac, 5e2 * dt_frac, label="dt")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("dt")
    ax.set_ylabel("Error in K")
    ax.set_ylim([1e-4, 1e2])
    ax.set_xlim([1e-4, 1e-1])
    plt.legend()
    plt.show()
    return err_temps


def show_Tprofile(cn_temps, im_temps, it, step):
    """
    Show the different temperature for each method at one step,
    one iteration.
    """
    nx = 101
    length = 2
    spaces = np.linspace(0, length, nx)
    spaces[1] = spaces[2] / 3
    # Temp
    cn_temp = cn_temps[step][:, it]
    im_temp = im_temps[step][:, it]
    ref_temp = cn_temps[4096][:, it * int(4096 / step)]
    # Plot
    fig, ax = plt.subplots()
    ax.plot(spaces, cn_temp, label="CN")
    ax.plot(spaces, im_temp, label="IM")
    ax.plot(spaces, ref_temp, label="Ref")
    ax.set_xscale("log")
    plt.legend()
    plt.show()


cn_times = np.load("cn_times.npy", allow_pickle=True).item()
cn_temps = np.load("cn_temps.npy", allow_pickle=True).item()
im_times = np.load("im_times.npy", allow_pickle=True).item()
im_temps = np.load("im_temps.npy", allow_pickle=True).item()

# Get params
steps = np.array(list(im_times.keys()))
ns = steps.size
dts = np.zeros(ns)
for isp, val in enumerate(im_times.values()):
    dts[isp] = val[0]

period = 86400 * 79.3
dt_frac = dts / period

im_errs = err_with_ref(im_temps, im_times, cn_temps[4096], cn_times[4096])
im_err_mean = get_err_mean(im_errs)
im_err_max = get_err_max(im_errs)
# Cn
cn_errs = err_with_ref(cn_temps, cn_times, cn_temps[4096], cn_times[4096])
cn_err_mean = get_err_mean(cn_errs)
cn_err_max = get_err_max(cn_errs)

fig, ax = plt.subplots()
ax.plot(dt_frac[1:], im_err_mean, "-o", label="IM Mean Error")
ax.plot(dt_frac[1:], cn_err_mean, "-o", label="CN Mean Error")
# ax.plot(dt_frac[1:], im_err_max, "-o", label="IM Max Error")
# ax.plot(dt_frac[1:], cn_err_max, "-o", label="CN Max Error")
ax.plot(dt_frac, 4e1 * dt_frac**2, "--", label="dt^2")
ax.plot(dt_frac, 1e1 * dt_frac, "--", label="dt")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$\Delta t/P$")
ax.set_ylabel("Error in K")
# ax.set_ylim([1e-4, 1e2])
# ax.set_xlim([1e-4, 1e-1])
plt.legend()
plt.savefig("conv.png")
plt.show()


# # it = 10
# # step = 32
# # show_Tprofile(cn_temps, im_temps, it, step)
