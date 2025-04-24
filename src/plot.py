import numpy as np
import matplotlib.pyplot as plt


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
    x_uniform = np.linspace(X.min(), X.max(), X.shape[1])
    y_uniform = np.linspace(Y.min(), Y.max(), Y.shape[0])
    X_uniform, Y_uniform = np.meshgrid(x_uniform, y_uniform)
    u_uniform = np.zeros_like(X_uniform)
    v_uniform = np.zeros_like(Y_uniform)

    for i in range(X.shape[0]):
        u_uniform[i, :] = np.interp(x_uniform, X[i, :], u[i, :])
    for j in range(X.shape[1]):
        v_uniform[:, j] = np.interp(y_uniform, Y[:, j], v[:, j])

    plt.streamplot(X_uniform, Y_uniform, u_uniform, v_uniform, density=3)
    plt.title('Flow Field')
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


def check_continuity(X, Y, u, v):
    dx = X[3, 6] - X[2, 6]
    dy = Y[6, 3] - Y[6, 2]

    du_dx = np.zeros_like(u)
    dv_dy = np.zeros_like(v)

    du_dx[:, 1:-1] = (u[:, 2:] - u[:, :-2]) / (2 * dx)
    dv_dy[1:-1, :] = (v[2:, :] - v[:-2, :]) / (2 * dy)

    divergence = du_dx + dv_dy

    max_divergence = np.max(np.abs(divergence[1:-1, 1:-1]))
    mean_divergence = np.mean(np.abs(divergence[1:-1, 1:-1]))
    print(f"Maximum divergence: {max_divergence:.2e}")
    print(f"Mean divergence: {mean_divergence:.2e}")


if __name__ == "__main__":
    # Load results from C++ simulation
    X, Y, psi, zeta, u, v = load_results("../buildDir/cavity_results.csv")

    # Plot results
    plot_results(Y, X, psi, zeta, u, v)
    check_continuity(X, Y, u, v)
