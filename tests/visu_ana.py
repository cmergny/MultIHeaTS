import matplotlib.pyplot as plt


def plot_compare(spaces, val_temps, sol_temps):
    fig, ax = plt.subplots()
    nt = val_temps.shape[1]
    n_curves = 5
    for it in range(0, nt, nt // n_curves):
        ax.plot(spaces, val_temps[:, it], color="blue")
        ax.plot(spaces, sol_temps[:, it], ".", color="red")
    # ax.set_ylim(0, 1)
    ax.set_xlabel("Depth (m)")
    ax.set_ylabel("Temperature (K)")
    plt.legend(["Analytic", "Solver"])
    plt.show()


def single_plot(spaces, temps):
    fig, ax = plt.subplots()
    for it in range(temps.shape[1]):
        ax.plot(spaces, temps[:, it], label=f"{it}")
    # ax.set_ylim(0, 1)
    ax.set_xlabel("Depth (m)")
    ax.set_ylabel("Temperature (K)")
    # plt.legend()
    plt.show()


def error_plot(dt, nt, err):
    fig, ax = plt.subplots()
    ax.loglog([dt * i for i in range(nt)], err, ".")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Sum of Errors (K)")
    plt.show()
