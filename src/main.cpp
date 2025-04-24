#include <omp.h>
#include <iostream>
#include "LidDrivenCavity.h"
#include "subprojects/tomlplusplus/toml.hpp"

int main(const int argc, char **argv) {
    omp_set_num_threads(6);
    std::string config_file = "config.toml";

    if (argc > 1) {
        config_file = argv[1];
    }
    try {
        auto config = toml::parse_file(config_file);

        const double Re = config["LidDrivenCavity"]["Re"].value_or(500.0);
        const double Length = config["LidDrivenCavity"]["length"].value_or(2.0);
        const int grid_points = config["LidDrivenCavity"]["grid_points"].value_or(100);
        const double time_step = config["LidDrivenCavity"]["time_step"].value_or(0.01);
        const int max_iterations = config["LidDrivenCavity"]["max_iterations"].value_or(1000);
        const double tolerance = config["LidDrivenCavity"]["tolerance"].value_or(1e-6);
        const bool python_plot = config["LidDrivenCavity"]["python_plot"].value_or(true);
        const double lid_velocity = config["LidDrivenCavity"]["lid_velocity"].value_or(1.0);

        LidDrivenCavity cavity(Re, Length, grid_points, time_step, max_iterations, tolerance, lid_velocity);
        cavity.solve();

        if (python_plot) {
            system("python plot.py");
        }
    } catch (const toml::parse_error &err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
        return 1;
    }

    return 0;
}
