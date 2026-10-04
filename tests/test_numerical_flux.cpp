#include "euler/eos.hpp"
#include "euler/numerical_flux.hpp"
#include "euler/state.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

// Consistency: with no jump across the face, the correction term vanishes and the
// Rusanov flux must reduce to the physical flux. A moving state (u != 0) exercises every term.
TEST_CASE("Rusanov flux of two identical states equals the physical flux", "[flux]") {
    const euler::Conserved q{2.0, 6.0, 19.0}; // (rho, u, p) = (2, 3, 4)
    const euler::IdealGas air{};

    const euler::Flux f_rusanov = euler::rusanov_flux(q, q, air);
    const euler::Flux f_physical = euler::physical_flux(euler::to_primitive(q, air), air);

    REQUIRE_THAT(f_rusanov.mass, Catch::Matchers::WithinRel(f_physical.mass, 1e-12));
    REQUIRE_THAT(f_rusanov.momentum, Catch::Matchers::WithinRel(f_physical.momentum, 1e-12));
    REQUIRE_THAT(f_rusanov.energy, Catch::Matchers::WithinRel(f_physical.energy, 1e-12));
}

// Sod diaphragm face at t = 0. Both sides are at rest, so a = max(c_L, c_R) = sqrt(1.4),
// the averaged mass and energy fluxes are 0, and the averaged momentum flux is (1 + 0.1) / 2.
TEST_CASE("Rusanov flux at the Sod diaphragm has the expected values", "[flux]") {
    const euler::Conserved left_u{1.0, 0.0, 2.5};
    const euler::Conserved right_u{0.125, 0.0, 0.25};
    const euler::IdealGas air{};

    const euler::Flux f = euler::rusanov_flux(left_u, right_u, air);

    REQUIRE_THAT(f.mass, Catch::Matchers::WithinRel(0.5 * std::sqrt(1.4) * 0.875, 1e-12));
    REQUIRE_THAT(f.momentum, Catch::Matchers::WithinRel(0.55, 1e-12));
    REQUIRE_THAT(f.energy, Catch::Matchers::WithinRel(0.5 * std::sqrt(1.4) * 2.25, 1e-12));
}
