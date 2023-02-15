import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


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


def plot_temp(spaces, temps, it):
    """Plot one temperature profile at one time."""
    fig, ax = plt.subplots()
    ax.plot(spaces, temps[it], label=f"{it}")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("Temperature (K)")
    ax.set_ylim(temps.min(), temps.max())
    ax.set_xscale("log")
    plt.legend()
    plt.show()


def animate_function(spaces, temps, step=1, save=False):
    """Plot an animation of the temperature with time"""
    fig, ax = plt.subplots(figsize=(5, 3.5))
    (line,) = ax.plot([])

    def animate(it):
        line.set_data((spaces, temps[it * step]))
        ax.set_title(f"{it}")
        return line

    ax.set_ylabel("Temperature (K)")
    ax.set_xlim(spaces.min(), spaces.max())
    ax.set_ylim(temps.min(), temps.max())
    ax.set_xscale("symlog")
    anim = FuncAnimation(fig, animate, frames=temps.shape[0], interval=1, repeat=True)
    return anim
