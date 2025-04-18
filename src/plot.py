import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def load_results(filename):
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)

    x = data[:, 0]
    y = data[:, 1]
    psi = data[:, 2]
    zeta = data[:, 3]
    u = data[:, 4]
    v = data[:, 5]
    N = int(np.sqrt(len(x)))
    X = x.reshape(N, N)
    Y = y.reshape(N, N)
    psi = psi.reshape(N, N)
    zeta = zeta.reshape(N, N)
    u = u.reshape(N, N)
    v = v.reshape(N, N)

    return X, Y, psi, zeta, u, v


def plot_results(X, Y, psi, zeta, u, v):
    plt.figure(figsize=(18, 6))

    # Stream function contour plot
    plt.subplot(1, 3, 1)
    plt.contourf(X, Y, psi, cmap='viridis')
    plt.colorbar(label='Stream Function')
    plt.title('Stream Function')
    plt.xlabel('x')
    plt.ylabel('y')

    # Vorticity contour plot
    plt.subplot(1, 3, 2)
    plt.contourf(X, Y, zeta, levels=20, cmap='coolwarm')
    plt.colorbar(label='Vorticity')
    plt.title('Vorticity')
    plt.xlabel('x')
    plt.ylabel('y')

    # Velocity quiver plot (subsampled for clarity)
    subsample = 5
    plt.subplot(1, 3, 3)
    plt.quiver(X[::subsample, ::subsample],
               Y[::subsample, ::subsample],
               u[::subsample, ::subsample],
               v[::subsample, ::subsample])
    plt.title('Velocity Field')
    plt.xlabel('x')
    plt.ylabel('y')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Load results from C++ simulation
    X, Y, psi, zeta, u, v = load_results("../buildDir/cavity_results.csv")

    # Plot results
    plot_results(Y, X, psi, zeta, u, v)
