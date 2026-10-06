// Runs the Sod shock tube with the first-order finite-volume solver and writes the result
// at t = 0.2 to sod_first_order.csv (columns: x, rho, u, p).

#include "euler/eos.hpp"
#include "euler/finite_volume.hpp"
#include "euler/state.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

int main() {
    // Problem: Sod shock tube on [0, 1], diaphragm at x = 0.5, air (gamma = 1.4).
    const euler::Grid grid{100, 0.0, 1.0, euler::Reconstruction::first_order};
    const euler::IdealGas air{};
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const double x_interface = 0.5;
    const double t_end = 0.2;
    const double cfl = 0.5;

    std::vector<euler::Conserved> u =
        euler::riemann_initial_state(grid, left, right, x_interface, air);
    std::vector<euler::Flux> face_flux(grid.n_cells + 1); // allocated once, reused every step

    // Time loop: a fresh CFL step every iteration, clipped so the last step lands on t_end.
    double t = 0.0;
    std::size_t steps = 0;
    while (t < t_end) {
        const double dt = std::min(euler::cfl_time_step(u, grid, cfl, air), t_end - t);
        euler::advance(u, face_flux, grid, dt, air);
        t += dt;
        ++steps;
    }

    // Output: one line per interior cell with its centre and primitive variables.
    const std::size_t g = euler::ghost_cells_per_side(grid.reconstruction);
    const double dx = euler::cell_width(grid);
    std::ofstream file("sod_first_order.csv");
    file << "x,rho,u,p\n";
    for (std::size_t i = 0; i < grid.n_cells; ++i) {
        const double x = grid.x_min + (static_cast<double>(i) + 0.5) * dx;
        const euler::Primitive w = euler::to_primitive(u[g + i], air);
        file << x << ',' << w.rho << ',' << w.u << ',' << w.p << '\n';
    }

    std::cout << "Sod: " << steps << " steps, t = " << t << '\n';
    return EXIT_SUCCESS;
}
