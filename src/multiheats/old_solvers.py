import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import shape


class HeatSolver:
    """A class object for numerical solving of the heat equation."""

    def __init__(self, nx, nt, rho, cp, k, dt, dx, L=1, ux_0=None) -> None:
        """Initialise the solver with heat properties
        and simulation parameters"""
        # Finite diff params
        self.nx = nx
        self.nt = nt
        self.L = L
        # Space/Time increments
        self.dx = dx
        self.dt = dt
        # Heat parameters
        self.rho = rho
        self.cp = cp
        self.k = k
        self.alpha = self.k / self.rho / self.cp
        self.r = self.dt / self.dx**2 / self.rho / self.cp
        # Temperature array
        self.u = np.zeros((self.nx, self.nt))
        self.u = np.zeros((self.nx, self.nt))
        # Set Initial Condition
        self.ux_0 = np.zeros(self.nx) if ux_0 is None else ux_0
        self.u[:, 0] = self.ux_0
        self.u[:, 0] = self.ux_0
        # Set Source term
        t = np.linspace(0, 10, nt)
        self.Q = np.zeros((self.nx, self.nt))
        # self.Q = np.random.uniform(size=(nx, nt))
        # Set BC
        self.BC = np.zeros((2, self.nt))
        # self.BC[:, :] = 1

    def plot(self, array=None):
        array = self.u if array is None else array
        """Plots 10 time iterations of the solution."""
        times = [k for k in range(0, self.nt, self.nt // 10)]
        fig, ax = plt.subplots()
        for t in times:
            ax.plot(
                np.linspace(0, self.L, self.nx),
                array[:, t],
                label=f"$t = {t}$",
                color="blue",
                alpha=0.3,
            )
        ax.set_xlabel("x")
        ax.set_ylabel(r"$u(x, t_0)$")
        # plt.legend([r'i_$th$ solution'])
        plt.show()

    def animate(self, array=None, save="False", frame_step=None):
        """Animate the evolution of u"""
        # Choose implicit sol if no array provided
        array = self.u if array is None else array
        frame_step = 10 if frame_step is None else frame_step
        # Define fct to call each frame
        def update(i):
            """Updates each frame"""
            for j in range(width):
                ar[:, j] = array[:, i * frame_step]
            ax.pcolor(ar, vmin=min_, vmax=max_, cmap="coolwarm")
            ax.set_title(f"frame {i}")
            ax.set_xlabel("Width")
            ax.set_ylabel("x")

        # Set vars
        width = 20
        ar = np.zeros((array.shape[0], width))
        min_, max_ = array.min(), array.max()
        # Animated plot
        fig, ax = plt.subplots()
        for j in range(width):
            ar[:, j] = array[:, 0]
        im = ax.pcolor(ar, vmin=min_, vmax=max_, cmap="coolwarm")
        ani = FuncAnimation(
            fig, update, frames=array.shape[1] // frame_step, interval=10, repeat=False
        )
        fig.colorbar(im).set_label("u(x)")
        if save:
            ani.save(
                "figures/heat_loop.gif",
                writer="imagemagick",
                progress_callback=lambda i, n: print(
                    "Progress: " + str(round(i / n * 100, 1)) + "%"
                ),
            )
        else:
            plt.show()
        return ani

    def stability_analysis(self):
        stab = (self.alpha * self.dt / self.dx**2).max()
        print(f"Stability number = {stab}")


class ExplicitSolver(HeatSolver):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.name = "ExplicitSolver"

    def explicit_scheme(self):
        """
        Solves the discretized heat equation explicitely with Euler Forward Scheme.
        Returns: u[:, :] np.array with u.shape = (nx, nt)
        """

        # Loop through space and time
        for it in range(1, self.nt):
            for ix in range(self.nx):
                r = (
                    self.dt
                    / self.dx**2
                    / (self.rho[ix, it - 1] * self.cp[ix, it - 1])
                )
                # Preliminary iterms
                term1, term2 = self.set_BC(ix, it, "Neumann")
                # Final terms
                self.u[ix, it] = self.u[ix, it - 1] + self.r[ix, it - 1] * (
                    term1 + term2 + self.Q[ix, it] * self.dx**2
                )

        return self.u

    def set_BC(self, ix, it, mode="Neumann"):
        N = self.nx - 1
        k = self.k
        if mode == "Neumann":
            if ix == 0:
                term1 = (
                    (self.k[1, it - 1] - self.k[0, it - 1])
                    * self.dx
                    * self.BC[0, it - 1]
                )
                term2 = (
                    2
                    * self.k[0, it - 1]
                    * (
                        self.u[1, it - 1]
                        - self.u[0, it - 1]
                        - self.dx * self.BC[0, it - 1]
                    )
                )
            elif ix == N:
                term1 = (
                    (self.k[N, it - 1] - self.k[N - 1, it - 1])
                    * self.dx
                    * self.BC[1, it - 1]
                )
                term2 = (
                    2
                    * self.k[N, it - 1]
                    * (
                        self.u[N - 1, it - 1]
                        - self.u[N, it - 1]
                        + self.dx * self.BC[1, it - 1]
                    )
                )
            else:
                term1 = (
                    (k[ix + 1, it - 1] - k[ix - 1, it - 1])
                    * (self.u[ix + 1, it - 1] - self.u[ix - 1, it - 1])
                    / 4
                )  # Divide by 4 important
                term2 = k[ix, it - 1] * (
                    self.u[ix + 1, it - 1]
                    - 2 * self.u[ix, it - 1]
                    + self.u[ix - 1, it - 1]
                )

        # Not working for now
        elif mode == "Dirichlet":
            if ix == 0:
                self.u[0, it] = self.BC[0, it]
            if ix == N:
                self.u[N, it] = self.BC[N, it]
        return (term1, term2)


class OldImplicitSolver(HeatSolver):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.name = "ImplicitSolver"

    def implicit_scheme(self):
        """
        Solves the discretized heat equation implicitely with Euler Backward Scheme.
        Returns: u[:, :] np.array with u.shape=(nx, nt)
        """
        r = self.r
        # Define matrix
        A = np.zeros((self.nx, self.nx, self.nt))
        # Loop through time
        for it in range(1, self.nt):
            for ix in range(1, self.nx - 1):
                A[ix, ix - 1] = (
                    r[ix, it - 1]
                    / 4
                    * (
                        self.k[ix + 1, it - 1]
                        - self.k[ix - 1, it - 1]
                        - 4 * self.k[ix, it - 1]
                    )
                )  # an
                A[ix, ix] = 1 + 2 * r[ix, it - 1] * self.k[ix, it - 1]  # bn
                A[ix, ix + 1] = (
                    r[ix, it - 1]
                    / 4
                    * (
                        self.k[ix - 1, it - 1]
                        - self.k[ix + 1, it - 1]
                        - 4 * self.k[ix, it - 1]
                    )
                )  # cn
            # Add source term
            source = (
                np.copy(self.u[:, it - 1])
                + r[ix, it - 1] * self.dx**2 * self.Q[ix, it - 1]
            )
            # Set BC
            A, source = self.set_BC(A, source, it)
            # Inverse matrix
            Ainv = np.linalg.inv(A[:, :, it])
            # Compute next iteration
            self.u[:, it] = np.dot(Ainv, source)
        return self.u

    def set_BC(self, A, source, it, mode="Neumman"):
        """Set boundary conditions for implicit Euler Scheme"""
        if mode == "Neumman":
            source[0] = self.u[0, it - 1] + self.r[0, it - 1] * (
                self.dx * self.BC[0, it] * (self.k[1, it - 1] - 3 * self.k[0, it - 1])
                + self.dx**2 * self.Q[0, it - 1]
            )
            source[-1] = self.u[-1, it - 1] + self.r[-1, it] * (
                self.dx
                * self.BC[-1, it]
                * (3 * self.k[-1, it - 1] - self.k[-2, it - 1])
                + self.dx**2 * self.Q[-1, it - 1]
            )
            # Source = T-1 si BC=0 et Q=0
            # Matrix
            A[0, 0, :] = 1 + 2 * self.r[0, it - 1] * self.k[0, it - 1]  # b1
            A[0, 1, :] = -2 * self.r[0, it - 1] * self.k[0, it - 1]  # c1
            A[self.nx - 1, self.nx - 2, :] = (
                -2 * self.r[-1, it - 1] * self.k[-1, it - 1]
            )  # an
            A[self.nx - 1, self.nx - 1, :] = (
                1 + 2 * self.r[-1, it - 1] * self.k[-1, it - 1]
            )  # bn
        # Change mode
        elif mode == "Dirichlet":
            source[0] = self.BC[0, it]
            source[-1] = self.BC[1, it]
            # Matrix
            A[0, 0, :] = 1
            A[0, 1, :] = 0
            A[self.nx - 1, self.nx - 2, :] = 0
            A[self.nx - 1, self.nx - 1, :] = 1
        return (A, source)


class OldValidator(HeatSolver):
    """Class used to validate the numerical solver solutions.
    It computes analytical solutions for known initial temperatures."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.name = "AnalyticSolution"
        self.alpha = self.alpha.mean()

    def solve_stepfunction(self):
        """Computes the analytic sol of the heat equa in [0,L] with:
             u(x, 0) = 0 if x<L/2 and 1 if x>L/2
             u'(O) = u'(L) = 0
        Returns: u: array of shape = (nx, nt)
        """
        self.name = "Step Function"
        # Fourier parameters
        N = 50  # nbr of fourier terms to compute
        a0 = 1  # first term
        un = np.zeros(N)
        self.u[self.nx // 2 :, :] = 1
        # Space/Time loop
        for it in range(1, self.nt):
            for ix in range(self.nx):
                x = self.dx * ix
                t = self.dt * it
                # Compute Fourier coefficients
                for n in range(1, N, 2):
                    an = 2 / np.pi / n * (-1) ** ((n - 1) / 2 + 1)
                    bn = 0
                    wn = n * np.pi * x / self.L  # pulsation
                    un[n] = (an * np.cos(wn) + bn * np.sin(wn)) * np.exp(
                        -self.alpha * (n * np.pi / self.L) ** 2 * t
                    )
                # Sum coeff for solution
                self.u[ix, it] = a0 / 2 + un.sum()
        return self.u
