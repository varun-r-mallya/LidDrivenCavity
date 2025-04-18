//
// Created by xeon on 4/19/25.
//

#ifndef LID_DRIVEN_CAVITY_H
#define LID_DRIVEN_CAVITY_H

#include <vector>
#include <chrono>


class LidDrivenCavity {
public:
    explicit LidDrivenCavity(double Re = 500.0, double L = 2.0, int N = 100,
                             double dt = 0.01, int max_iter = 10000, double tol = 1e-6);

    void solve();

    void save_results(const std::string &filename) const;

private:
    void apply_boundary_conditions();

    void solve_vorticity();

    void solve_stream_function();

    void compute_velocities();

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
    std::vector<std::vector<double> > vorticity; // Vorticity
    std::vector<std::vector<double> > u; // x-velocity
    std::vector<std::vector<double> > v; // y-velocity

    // Grid to set stuff inside the CSV
    std::vector<double> x;
    std::vector<double> y;
};


#endif //LID_DRIVEN_CAVITY_H
