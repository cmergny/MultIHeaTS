import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def animate_function(spaces, temps, interf=0, step=1, frames=None, save=False):
    """Plot an animation of the temperature with time"""
    frames = temps.shape[0] // step if frames is None else frames
    fig, ax = plt.subplots()
    (line,) = ax.plot([])
    ax.axvline(x=interf, alpha=0.6, linestyle="--", color="grey", label="Interface")

    def animate(it):
        line.set_data((spaces, temps[it * step]))
        ax.set_title(f"{it}")
        return line

    ax.set_ylabel("Temperature (K)")
    ax.set_xlabel("Depth (m)")
    ax.set_ylim(temps.min(), temps.max())
    ax.set_xlim([spaces[1], spaces.max()])
    ax.set_xscale("log")
    anim = FuncAnimation(fig, animate, frames=frames, interval=50, repeat=False)
    plt.title("Temperature Evolution")
    plt.legend()

    if save:
        progress_callback = lambda i, n: print(f"Saving frame {i} of {n}")
        anim.save("../../figures/temp_evo.mp4", progress_callback=progress_callback)
    return anim


def animate_2functions(
    spaces,
    a_temps,
    b_temps,
    a_label="A",
    b_label="B",
    interf=0,
    step=1,
    frames=None,
    save=False,
):
    """Plot an animation of the temperature with time"""
    frames = a_temps.shape[0] // step if frames is None else frames
    fig, ax = plt.subplots()
    (aline,) = ax.plot([], color="red", label=a_label)
    (bline,) = ax.plot([], color="blue", label=b_label)
    ax.axvline(x=interf, alpha=0.6, linestyle="--", color="grey", label="Interface")

    def animate(it):
        aline.set_data((spaces, a_temps[it * step]))
        bline.set_data((spaces, b_temps[it * step]))
        ax.set_title(f"{it}")
        return aline, bline

    ax.set_ylabel("Temperature (K)")
    ax.set_xlabel("Depth (m)")
    ax.set_ylim(a_temps.min(), a_temps.max())
    ax.set_xlim([spaces[1], spaces.max()])
    ax.set_xscale("log")
    anim = FuncAnimation(fig, animate, frames=frames, interval=50, repeat=False)
    plt.title("Temperature Evolution")
    plt.legend()

    if save:
        progress_callback = lambda i, n: print(f"Saving frame {i} of {n}")
        anim.save("temp_evo.mp4", progress_callback=progress_callback)
    return anim
