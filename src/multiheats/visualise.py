import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def use_latex():
    tex_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 13,
        "font.size": 11,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
    }
    plt.rcParams.update(tex_fonts)


def plot_multi_temp(spaces, temps, n_curves=10):
    """Plot multiple temperature profiles chosen evenly among all times."""
    fig, ax = plt.subplots()
    nt = temps.shape[0]
    for it in range(0, nt, nt // n_curves):
        ax.plot(spaces, temps[it], label=f"{it}")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("Temperature (K)")
    # ax.set_xscale("log")
    plt.legend()
    plt.show()


def plot_temp(spaces, temps, it, interf=0):
    """Plot one temperature profile at one time."""
    mean_temp = temps[it].mean()
    fig, ax = plt.subplots()
    ax.axvline(x=interf, alpha=0.2, linestyle="-", color="black", label="Interface")
    ax.plot(spaces, temps[it])
    ax.set_xlabel("Depth (m)")
    ax.set_ylabel("Temperature (K)")
    ax.set_ylim(0.9 * mean_temp, 1.1 * mean_temp)
    ax.set_xlim(0, spaces.max())
    ax.set_xscale("symlog")
    plt.title("Bilayer temperature profile snapshot")
    plt.legend()
    # plt.show()
    plt.savefig("../../examples/temp_bilayer.png")


def animate_function(spaces, temps, interf=0, step=1, frames=None, save=False):
    """Plot an animation of the temperature with time"""
    frames = temps.shape[0] // step if frames is None else frames
    fig, ax = plt.subplots()
    (line,) = ax.plot([])
    ax.axvline(x=interf, alpha=0.3, linestyle="-", color="black")

    def animate(it):
        line.set_data((spaces, temps[it * step]))
        # ax.set_title(f"{it}")
        return line

    ax.set_ylabel("Temperature (K)")
    ax.set_xlabel("Depth (m)")
    ax.set_xlim(spaces.min(), spaces.max())
    ax.set_ylim(temps.min(), temps.max())
    ax.set_xscale("symlog")
    anim = FuncAnimation(fig, animate, frames=frames, interval=1, repeat=False)
    plt.title("Temperature Evolution")

    if save:
        progress_callback = lambda i, n: print(f"Saving frame {i} of {n}")
        anim.save("temp_evo.gif", progress_callback=progress_callback)
    return anim
