#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include "renderer.h"
#include <omp.h>

class LidDrivenCavity {
public:
    LidDrivenCavity(double Re = 100.0, double L = 1.0, int N = 100,
                    double dt = 0.001, int max_iter = 10000, double tol = 1e-6)
        : Re(Re), L(L), N(N), dt(dt), max_iter(max_iter), tol(tol),
          dx(L / (N - 1)), dy(L / (N - 1)), u_top(1.0) {

        psi = std::vector<std::vector<double> >(N, std::vector<double>(N, 0.0));
        zeta = std::vector<std::vector<double> >(N, std::vector<double>(N, 0.0));
        u = std::vector<std::vector<double> >(N, std::vector<double>(N, 0.0));
        v = std::vector<std::vector<double> >(N, std::vector<double>(N, 0.0));

        // Initialize grid
        x = std::vector<double>(N);
        y = std::vector<double>(N);
#pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            x[i] = i * dx;
            y[i] = i * dy;
        }
    }

    void solve() {
        auto start = std::chrono::high_resolution_clock::now();

        Renderer renderer(800, 800, u, v, x, y);
        renderer.initialize();

        // #pragma omp parallel for reduction(max:max_diff)
        for (int n = 0; n < max_iter; ++n) {
            // Store old vorticity for convergence check
            auto zeta_old = zeta;
            // Check if window should close
            renderer.handleEvents();
            if (!renderer.isWindowOpen()) {
                break;
            }

            // Perform one iteration
            apply_boundary_conditions();
            solve_vorticity();
            solve_streamfunction();
            compute_velocities();

            renderer.updateData(u, v);
            renderer.render();

            if (n % 100 == 0) {
                std::cout << "Iteration: " << n << std::endl;
            }
            // Check for convergence
            double max_diff = 0.0;
#pragma omp parallel for reduction(+:max_diff)
            for (int i = 0; i < N; ++i) {
#pragma omp parallel for
                for (int j = 0; j < N; ++j) {
                    double diff = std::abs(zeta[i][j] - zeta_old[i][j]);
                    if (diff > max_diff) max_diff = diff;
                }
            }

            if (n % 100 == 0) {
                std::cout << "Iteration: " << n << ", Max change: " << max_diff << std::endl;
                save_results("cavity_results.csv");
                // system("python plot.py");
            }

            if (max_diff < tol) {
                std::cout << "Converged after " << n << " iterations" << std::endl;
                break;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Simulation completed in " << duration.count() << " ms" << std::endl;

        while (renderer.isWindowOpen()) {
            renderer.handleEvents();
            renderer.render();
        }
    }

    void save_results(const std::string &filename) {
        std::ofstream outfile(filename);

        // Write header
        outfile << "x,y,psi,zeta,u,v\n";

        // Write data
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                outfile << x[i] << "," << y[j] << ","
                        << psi[i][j] << "," << zeta[i][j] << ","
                        << u[i][j] << "," << v[i][j] << "\n";
            }
        }

        outfile.close();
        std::cout << "Results saved to " << filename << std::endl;
    }

private:
    void apply_boundary_conditions() {
        // Top wall (moving lid)
#pragma omp parallel for
        for (int j = 1; j < N - 1; ++j) {
            zeta[N - 1][j] = -2.0 * psi[N - 2][j] / (dy * dy) - 2.0 * u_top / dy;
        }

        // Bottom wall
#pragma omp parallel for
        for (int j = 1; j < N - 1; ++j) {
            zeta[0][j] = -2.0 * psi[1][j] / (dy * dy);
        }

        // Left wall
#pragma omp parallel for
        for (int i = 1; i < N - 1; ++i) {
            zeta[i][0] = -2.0 * psi[i][1] / (dx * dx);
        }

        // Right wall
#pragma omp parallel for
        for (int i = 1; i < N - 1; ++i) {
            zeta[i][N - 1] = -2.0 * psi[i][N - 2] / (dx * dx);
        }

        // Corners (average of adjacent sides)
        zeta[0][0] = 0.5 * (zeta[0][1] + zeta[1][0]);
        zeta[0][N - 1] = 0.5 * (zeta[0][N - 2] + zeta[1][N - 1]);
        zeta[N - 1][0] = 0.5 * (zeta[N - 1][1] + zeta[N - 2][0]);
        zeta[N - 1][N - 1] = 0.5 * (zeta[N - 1][N - 2] + zeta[N - 2][N - 1]);
    }

    void solve_vorticity() {
        auto zeta_new = zeta; // Create copy
#pragma omp parallel for
        for (int i = 1; i < N - 1; ++i) {
#pragma omp parallel for
            for (int j = 1; j < N - 1; ++j) {
                // Convection terms
                double u_zeta_x = u[i][j] * (zeta[i][j + 1] - zeta[i][j - 1]) / (2.0 * dx);
                double v_zeta_y = v[i][j] * (zeta[i + 1][j] - zeta[i - 1][j]) / (2.0 * dy);

                // Diffusion terms
                double diff_x = (zeta[i][j + 1] - 2.0 * zeta[i][j] + zeta[i][j - 1]) / (dx * dx);
                double diff_y = (zeta[i + 1][j] - 2.0 * zeta[i][j] + zeta[i - 1][j]) / (dy * dy);

                // Update vorticity
                zeta_new[i][j] = zeta[i][j] + dt * (
                                     -u_zeta_x - v_zeta_y + (1.0 / Re) * (diff_x + diff_y)
                                 );
            }
        }

        zeta = zeta_new;
    }

    void solve_streamfunction() {
        const double omega = 1.5; // Relaxation parameter
#pragma omp parallel for
        for (int k = 0; k < 100; ++k) {
#pragma omp parallel for
            for (int i = 1; i < N - 1; ++i) {
#pragma omp parallel for
                for (int j = 1; j < N - 1; ++j) {
                    // Update stream function using SOR
                    double psi_new = (1.0 - omega) * psi[i][j] + omega * (
                                         (dy * dy * (psi[i][j + 1] + psi[i][j - 1]) +
                                          dx * dx * (psi[i + 1][j] + psi[i - 1][j]) +
                                          dx * dx * dy * dy * zeta[i][j]) /
                                         (2.0 * (dx * dx + dy * dy))
                                     );
                    psi[i][j] = psi_new;
                }
            }
        }
    }

    void compute_velocities() {
        // Interior points
#pragma omp parallel for
        for (int i = 1; i < N - 1; ++i) {
#pragma omp parallel for
            for (int j = 1; j < N - 1; ++j) {
                u[i][j] = (psi[i + 1][j] - psi[i - 1][j]) / (2.0 * dy);
                v[i][j] = -(psi[i][j + 1] - psi[i][j - 1]) / (2.0 * dx);
            }
        }

        // Boundary conditions for velocity
#pragma omp parallel for
        for (int j = 0; j < N; ++j) {
            u[0][j] = 0.0; // Bottom wall
            u[N - 1][j] = u_top; // Top wall (moving lid)
            v[0][j] = 0.0; // Bottom wall
            v[N - 1][j] = 0.0; // Top wall
        }
#pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            u[i][0] = 0.0; // Left wall
            u[i][N - 1] = 0.0; // Right wall
            v[i][0] = 0.0; // Left wall
            v[i][N - 1] = 0.0; // Right wall
        }
    }

    // Parameters
    double Re; // Reynolds number
    double L; // Domain size
    int N; // Grid points
    double dt; // Time step
    int max_iter; // Maximum iterations
    double tol; // Tolerance

    // Derived parameters
    double dx, dy; // Grid spacing
    double u_top; // Lid velocity

    // Arrays
    std::vector<std::vector<double> > psi; // Stream function
    std::vector<std::vector<double> > zeta; // Vorticity
    std::vector<std::vector<double> > u; // x-velocity
    std::vector<std::vector<double> > v; // y-velocity

    // Grid
    std::vector<double> x;
    std::vector<double> y;
};

int main() {
    // Create and run solve
    omp_set_num_threads(omp_get_max_threads());
    LidDrivenCavity cavity(1000.0, 1.0, 100, 0.01, 1000);

    cavity.solve();
    // cavity.save_results("cavity_results.csv");

    return 0;
}
