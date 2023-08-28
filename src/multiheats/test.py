fig, axs = plt.subplots(1, 2, figsize=(9, 4))

axs[0].plot(dt_fracs, err_maxs, "-", marker=".", label="Implicit", color="#AA4499")
# axs[0].plot(dt_fracs[0], spen_max, "x", label="Explicit", color="#44AA99")
axs[0].set_xlabel(r"$\Delta t/P$")
axs[0].set_ylabel("Max Error (K)")
axs[0].set_yscale("asinh")
axs[0].set_xscale("log")
axs[0].set_ylim([3e-2, 50])
axs[0].axvspan(0, 0.5, alpha=0.2, color="blue", label="Explicit Stability Zone")
# rect = matplotlib.patches.Rectangle((0,0), 99, 10, color="orange", alpha=0.2, label="Accurate Computation Zone")
# axs[0].add_patch(rect)
axs[0].legend()

ax1 = axs[0].twiny()
ax1.plot(frac, np.zeros(day_frac.size))
ax1.set_xscale("log")
ax1.set_xlabel("Timestep as fraction of surface flux period")
ax1.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(decimals=2))


axs[1].plot(dt_fracs, err_means, "-", marker=".", label="Implicit", color="#AA4499")
axs[1].plot(
    dt_fracs, cn_means, "--", marker=".", label="Crank-Nicolson", color="#999933"
)
axs[1].plot(dt_fracs[0], spen_mean, "x", label="Explicit", color="#44AA99")
axs[1].set_xlabel("Fourier Number $F$")
axs[1].set_ylabel("Mean Error (K)")
axs[1].set_yscale("log")
axs[1].set_xscale("log")
axs[1].axvspan(0, 0.5, alpha=0.2, color="blue", label="Explicit Stability Zone")
axs[1].set_ylim(3e-2, 2.5)

ax2 = axs[1].twiny()
ax2.plot(frac, np.zeros(day_frac.size))
ax2.set_xscale("log")
ax2.set_xlabel("Timestep as fraction of surface flux period")
ax2.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(decimals=2))

axs[1].legend()
