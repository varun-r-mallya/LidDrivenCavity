//
// Created by xeon on 4/19/25.
//

#include "LidDrivenCavity.h"

#include "renderer.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <fstream>


auto LidDrivenCavity::compute_velocities() -> void {
    // Interior points
#pragma omp parallel for
    for (int i = 1; i < N - 1; ++i) {
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

void LidDrivenCavity::solve_stream_function() {
    constexpr double relaxation_factor = 0.2;
#pragma omp parallel for
    for (int k = 0; k < 30; ++k) {
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                // Update stream function using SOR
                const double psi_new = (1.0 - relaxation_factor) * psi[i][j] + relaxation_factor * (
                                           (dy * dy * (psi[i][j + 1] + psi[i][j - 1]) +
                                            dx * dx * (psi[i + 1][j] + psi[i - 1][j]) +
                                            dx * dx * dy * dy * vorticity[i][j]) /
                                           (2.0 * (dx * dx + dy * dy))
                                       );
                psi[i][j] = psi_new;
            }
        }
    }
}

auto LidDrivenCavity::solve_vorticity() -> void {
    auto zeta_new = vorticity; // Create copy

    // this is an FTCS scheme with the following stability condition
    // const double nu = u_top * L / Re;
    // assert((nu*dt/(dx*dx) + nu*dt/(dy*dy)) > 0.5);
#pragma omp parallel for
    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            // Convection terms
            const double u_zeta_x = u[i][j] * (vorticity[i][j + 1] - vorticity[i][j - 1]) / (2.0 * dx);
            const double v_zeta_y = v[i][j] * (vorticity[i + 1][j] - vorticity[i - 1][j]) / (2.0 * dy);

            // Diffusion terms
            const double diff_x = (vorticity[i][j + 1] - 2.0 * vorticity[i][j] + vorticity[i][j - 1]) / (dx * dx);
            const double diff_y = (vorticity[i + 1][j] - 2.0 * vorticity[i][j] + vorticity[i - 1][j]) / (dy * dy);

            // Update vorticity
            zeta_new[i][j] = vorticity[i][j] + dt * (
                                 -u_zeta_x - v_zeta_y + (u_top * L / Re) * (diff_x + diff_y)
                             );
        }
    }

    vorticity = zeta_new;
}

auto LidDrivenCavity::apply_boundary_conditions() -> void {
    // First order boundary conditions
    // Top wall: ω_{i,N}
#pragma omp parallel for
    for (int j = 1; j < N - 1; ++j) {
        vorticity[N - 1][j] = 2.0 * (psi[N - 1][j] - psi[N - 2][j]) / (dy * dy) - 2.0 * u_top / dy;
    }

    // Bottom wall: ω_{i,1}
    constexpr double u_bottom = 0;
#pragma omp parallel for
    for (int j = 1; j < N - 1; ++j) {
        vorticity[0][j] = 2.0 * (psi[0][j] - psi[1][j]) / (dy * dy) + 2.0 * u_bottom / dy;
    }

    // Left wall: ω_{1,j}
    constexpr double v_left = 0;
#pragma omp parallel for
    for (int i = 1; i < N - 1; ++i) {
        vorticity[i][0] = 2.0 * (psi[i][0] - psi[i][1]) / (dx * dx) - 2.0 * v_left / dx;
    }

    // Right wall: ω_{M,j}
    constexpr double v_right = 0;
#pragma omp parallel for
    for (int i = 1; i < N - 1; ++i) {
        vorticity[i][N - 1] = 2.0 * (psi[i][N - 1] - psi[i][N - 2]) / (dx * dx) + 2.0 * v_right / dx;
    }

    //         //second order boundary conditions
    //         // Top wall: ω_{i,N}
    // #pragma omp parallel for
    //         for (int j = 1; j < N - 1; ++j) {
    //             vorticity[N - 1][j] = (-7.0 * psi[N - 1][j] + 8.0 * psi[N - 2][j] - psi[N - 3][j]) / (dy * dy)
    //                              - 3.0 * u_top / dy;
    //         }
    //
    //         // Bottom wall: ω_{i,1}
    // #pragma omp parallel for
    //         for (int j = 1; j < N - 1; ++j) {
    //             constexpr double u_bottom = 0;
    //             vorticity[0][j] = (7.0 * psi[0][j] - 8.0 * psi[1][j] + psi[2][j]) / (dy * dy)
    //                          + 3.0 * u_bottom / dy;
    //         }
    //
    //         // Left wall: ω_{1,j}
    // #pragma omp parallel for
    //         for (int i = 1; i < N - 1; ++i) {
    //             constexpr double v_left = 0;
    //             vorticity[i][0] = (7.0 * psi[i][0] - 8.0 * psi[i][1] + psi[i][2]) / (dx * dx)
    //                          - 3.0 * v_left / dx;
    //         }
    //
    //         // Right wall: ω_{M,j}
    // #pragma omp parallel for
    //         for (int i = 1; i < N - 1; ++i) {
    //             constexpr double v_right = 0;
    //             vorticity[i][N - 1] = (-7.0 * psi[i][N - 1] + 8.0 * psi[i][N - 2] - psi[i][N - 3]) / (dx * dx)
    //                              + 3.0 * v_right / dx;
    //         }


    // Corners (average of adjacent sides)
    vorticity[0][0] = 0.5 * (vorticity[0][1] + vorticity[1][0]);
    vorticity[0][N - 1] = 0.5 * (vorticity[0][N - 2] + vorticity[1][N - 1]);
    vorticity[N - 1][0] = 0.5 * (vorticity[N - 1][1] + vorticity[N - 2][0]);
    vorticity[N - 1][N - 1] = 0.5 * (vorticity[N - 1][N - 2] + vorticity[N - 2][N - 1]);
}

LidDrivenCavity::LidDrivenCavity(const double Re, const double L, const int N,
                                 const double dt, const int max_iter, const double tol)
    : Re(Re), L(L), N(N), dt(dt), max_iter(max_iter), tol(tol),
      dx(L / (N - 1)), dy(L / (N - 1)), u_top(1.0) {
    psi = std::vector<std::vector<double> >(N, std::vector<double>(N, 0.0));
    vorticity = std::vector<std::vector<double> >(N, std::vector<double>(N, 0.0));
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

auto LidDrivenCavity::solve() -> void {
    const auto start = std::chrono::high_resolution_clock::now();

    Renderer renderer(800, 800, u, v, x, y);
    renderer.initialize();

    for (int n = 0; n < max_iter; ++n) {
        // Store old vorticity for convergence check
        auto zeta_old = vorticity;
        // Check if window should close
        renderer.handleEvents();
        if (!renderer.isWindowOpen()) {
            break;
        }

        // Perform one iteration
        apply_boundary_conditions();
        solve_vorticity();
        solve_stream_function();
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
            for (auto j = 0; j < N; ++j) {
                const double diff = std::abs(vorticity[i][j] - zeta_old[i][j]);
                if (diff > max_diff) max_diff = diff;
            }
        }

        if (n % 100 == 0) {
            std::cout << "Iteration: " << n << ", Max change: " << max_diff << std::endl;
            save_results("cavity_results.csv");
        }

        if (max_diff < tol) {
            std::cout << "Converged after " << n << " iterations" << std::endl;
            break;
        }
    }

    const auto end = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Simulation completed in " << duration.count() << " ms" << std::endl;

    while (renderer.isWindowOpen()) {
        renderer.handleEvents();
        renderer.render();
    }
}

auto LidDrivenCavity::save_results(const std::string &filename) const -> void {
    std::ofstream outfile(filename);

    // Write header
    outfile << "x,y,psi,vorticity,u,v\n";

    // Write data
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            outfile << x[i] << "," << y[j] << ","
                    << psi[i][j] << "," << vorticity[i][j] << ","
                    << u[i][j] << "," << v[i][j] << "\n";
        }
    }

    outfile.close();
    std::cout << "Results saved to " << filename << "\n";
}
