#include "euler/eos.hpp"
#include "euler/finite_volume.hpp"
#include "euler/state.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
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

// 4 cells on [0, 1] with one ghost per side: centres -0.125 | 0.125 0.375 0.625 0.875 | 1.125,
// so positions 0-2 get the left state and positions 3-5 the right state.
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
