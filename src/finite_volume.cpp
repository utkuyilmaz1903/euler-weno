#include "euler/finite_volume.hpp"

#include "euler/eos.hpp"
#include "euler/state.hpp"

#include <cstddef>

namespace euler {

std::size_t ghost_cells_per_side(Reconstruction reconstruction) noexcept {
    switch (reconstruction) {
    case Reconstruction::first_order:
        return 1;
    case Reconstruction::weno5:
        return 3;
    }
    return 0; // not reached: every Reconstruction is handled above
}

double cell_width(Grid grid) noexcept {
    const double dx = (grid.x_max - grid.x_min) / grid.n_cells;
    return dx;
}

std::vector<Conserved> riemann_initial_state(Grid grid, Primitive left, Primitive right,
                                             double x_interface, IdealGas gas) {
    const std::size_t g = ghost_cells_per_side(grid.reconstruction);
    const double dx = cell_width(grid);
    const std::size_t size = grid.n_cells + 2 * g;

    std::vector<Conserved> u(size);

    const Conserved u_left = to_conserved(left, gas);
    const Conserved u_right = to_conserved(right, gas);

    // j is the position in the list; j - g is the interior cell number (negative for the
    // left ghosts), so the centre is computed in double to avoid unsigned wrap-around.
    for (std::size_t j = 0; j < size; ++j) {
        const double x = grid.x_min + (static_cast<double>(j) - static_cast<double>(g) + 0.5) * dx;

        if (x < x_interface) {
            u[j] = u_left;
        } else {
            u[j] = u_right;
        }
    }
    return u;
}

void apply_transmissive_boundaries(std::span<Conserved> u, Grid grid) noexcept {
    const std::size_t g = ghost_cells_per_side(grid.reconstruction);

    const std::size_t first = g;
    const std::size_t last = u.size() - (1 + g);

    for (std::size_t k = 0; k < g; ++k) {
        u[k] = u[first];
        u[u.size() - (1 + k)] = u[last];
    }
}

} // namespace euler
