#pragma once

#include "euler/eos.hpp"
#include "euler/state.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace euler {

// How the solver gets the values on both sides of a face from the cell values.
// It decides how many neighbouring cells a face needs, and so how many ghost cells
// the grid needs at each end of the tube.
enum class Reconstruction {
    first_order, // a face uses only its two adjacent cells (Step 3, Rusanov)
    weno5,       // a face uses a five-cell stencil (Step 4)
};

// Number of ghost cells needed at each end of the tube for a reconstruction.
[[nodiscard]] std::size_t ghost_cells_per_side(Reconstruction reconstruction) noexcept;

// A tube [x_min, x_max] split into n_cells equal cells.
struct Grid {
    std::size_t n_cells;
    double x_min;
    double x_max;
    Reconstruction reconstruction;
};

// Width of one cell: dx = (x_max - x_min) / n_cells.
[[nodiscard]] double cell_width(Grid grid) noexcept;

// Initial cell values of a Riemann problem: a cell gets `left` if its centre lies left of
// x_interface and `right` otherwise. The returned list includes the ghost cells:
//   [ ghost cells | n_cells interior cells | ghost cells ]
[[nodiscard]] std::vector<Conserved>
riemann_initial_state(Grid grid, Primitive left, Primitive right, double x_interface, IdealGas gas);

// Transmissive (open-end) boundaries: every ghost cell gets a copy of the nearest interior
// cell, so waves leave the tube without reflecting. Called before every time step, because
// the interior cells change from step to step. `u` is the full list, ghosts included.
void apply_transmissive_boundaries(std::span<Conserved> u, Grid grid) noexcept;

} // namespace euler
