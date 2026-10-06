#include "euler/eos.hpp"
#include "euler/finite_volume.hpp"
#include "euler/state.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <vector>

TEST_CASE("Each reconstruction asks for the right number of ghost cells", "[fv]") {
    const std::size_t first_order_ghost_cells =
        euler::ghost_cells_per_side(euler::Reconstruction::first_order);
    const std::size_t weno5_ghost_cells = euler::ghost_cells_per_side(euler::Reconstruction::weno5);
    REQUIRE(first_order_ghost_cells == 1);
    REQUIRE(weno5_ghost_cells == 3);
}

TEST_CASE("Cell width is the tube length divided by the number of cells", "[fv]") {
    const euler::Grid grid{4, 0.0, 1.0, euler::Reconstruction::weno5};
    const double dx = euler::cell_width(grid);
    REQUIRE_THAT(dx, Catch::Matchers::WithinRel(0.25, 1e-12));
}

TEST_CASE("Sod initial state fills every cell, ghosts included", "[fv]") {
    const euler::Grid grid{4, 0.0, 1.0, euler::Reconstruction::first_order};
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    const std::vector<euler::Conserved> u =
        euler::riemann_initial_state(grid, left, right, 0.5, air);

    REQUIRE(u.size() == 6);
    for (const std::size_t j : {0, 2}) { // left ghost and last left cell
        REQUIRE_THAT(u[j].rho, Catch::Matchers::WithinRel(1.0, 1e-12));
        REQUIRE_THAT(u[j].E, Catch::Matchers::WithinRel(2.5, 1e-12));
    }
    for (const std::size_t j : {3, 5}) { // first right cell and right ghost
        REQUIRE_THAT(u[j].rho, Catch::Matchers::WithinRel(0.125, 1e-12));
        REQUIRE_THAT(u[j].E, Catch::Matchers::WithinRel(0.25, 1e-12));
    }
}

TEST_CASE("A WENO5 grid gets three ghost cells per side", "[fv]") {
    const euler::Grid grid{4, 0.0, 1.0, euler::Reconstruction::weno5};
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    const std::vector<euler::Conserved> u =
        euler::riemann_initial_state(grid, left, right, 0.5, air);
    REQUIRE(u.size() == 10); // 4 interior cells + 2 * 3 ghosts
}

// The interior end cells are changed by hand, as a time step would; the ghosts must follow.
TEST_CASE("Transmissive boundaries copy the end cells into the ghosts (first order)", "[fv]") {
    const euler::Grid grid{4, 0.0, 1.0, euler::Reconstruction::first_order};
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    std::vector<euler::Conserved> u = euler::riemann_initial_state(grid, left, right, 0.5, air);

    u[1] = euler::Conserved{2.0, 0.0, 5.0};
    u[4] = euler::Conserved{0.5, 0.0, 1.0};
    euler::apply_transmissive_boundaries(u, grid);

    REQUIRE_THAT(u[0].rho, Catch::Matchers::WithinRel(2.0, 1e-12));
    REQUIRE_THAT(u[0].E, Catch::Matchers::WithinRel(5.0, 1e-12));
    REQUIRE_THAT(u[5].rho, Catch::Matchers::WithinRel(0.5, 1e-12));
    REQUIRE_THAT(u[5].E, Catch::Matchers::WithinRel(1.0, 1e-12));
}

TEST_CASE("Transmissive boundaries fill all three ghost layers (WENO5)", "[fv]") {
    const euler::Grid grid{4, 0.0, 1.0, euler::Reconstruction::weno5};
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    std::vector<euler::Conserved> u = euler::riemann_initial_state(grid, left, right, 0.5, air);

    u[3] = euler::Conserved{2.0, 0.0, 5.0};
    u[6] = euler::Conserved{0.5, 0.0, 1.0};
    euler::apply_transmissive_boundaries(u, grid);

    for (std::size_t k = 0; k < 3; ++k) {
        REQUIRE_THAT(u[k].rho, Catch::Matchers::WithinRel(2.0, 1e-12));
        REQUIRE_THAT(u[k].E, Catch::Matchers::WithinRel(5.0, 1e-12));
        REQUIRE_THAT(u[9 - k].rho, Catch::Matchers::WithinRel(0.5, 1e-12));
        REQUIRE_THAT(u[9 - k].E, Catch::Matchers::WithinRel(1.0, 1e-12));
    }
}

// Sod at rest: the fastest wave is the sound speed on the left, c = sqrt(1.4).
TEST_CASE("CFL time step uses the fastest wave speed in the tube", "[fv]") {
    const euler::Grid grid{4, 0.0, 1.0, euler::Reconstruction::first_order};
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    const std::vector<euler::Conserved> u =
        euler::riemann_initial_state(grid, left, right, 0.5, air);

    const double dt = euler::cfl_time_step(u, grid, 0.5, air);
    REQUIRE_THAT(dt, Catch::Matchers::WithinRel(0.5 * 0.25 / std::sqrt(1.4), 1e-12));
}

// The fastest cell is the last interior cell, so the loop must include it.
TEST_CASE("CFL time step looks at every interior cell, including the last one", "[fv]") {
    const euler::Grid grid{4, 0.0, 1.0, euler::Reconstruction::first_order};
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    std::vector<euler::Conserved> u = euler::riemann_initial_state(grid, left, right, 0.5, air);

    // Layout: ghost 0 | interior 1 2 3 4 | ghost 5. Cell 4 moves at u = 2, so |u| + c = 2 +
    // sqrt(1.4).
    u[4] = euler::to_conserved(euler::Primitive{1.0, 2.0, 1.0}, air);

    const double dt = euler::cfl_time_step(u, grid, 0.5, air);
    REQUIRE_THAT(dt, Catch::Matchers::WithinRel(0.5 * 0.25 / (2.0 + std::sqrt(1.4)), 1e-12));
}
