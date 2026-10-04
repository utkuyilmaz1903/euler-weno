#include "euler/eos.hpp"
#include "euler/state.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

TEST_CASE("Sod left state converts to the expected conserved variables", "[eos]") {
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::IdealGas air{};
    const euler::Conserved left_u = euler::to_conserved(left, air);
    REQUIRE_THAT(left_u.rho, Catch::Matchers::WithinRel(1.0, 1e-12));
    REQUIRE_THAT(left_u.mom, Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(left_u.E, Catch::Matchers::WithinRel(2.5, 1e-12));
}

TEST_CASE("Sod right state converts to the expected conserved variables", "[eos]") {
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    const euler::Conserved right_u = euler::to_conserved(right, air);
    REQUIRE_THAT(right_u.rho, Catch::Matchers::WithinRel(0.125, 1e-12));
    REQUIRE_THAT(right_u.mom, Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(right_u.E, Catch::Matchers::WithinRel(0.25, 1e-12));
}

TEST_CASE("Sod left state has the expected sound speed", "[eos]") {
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::IdealGas air{};
    REQUIRE_THAT(euler::sound_speed(left, air), Catch::Matchers::WithinRel(std::sqrt(1.4), 1e-12));
}

TEST_CASE("Sod right state has the expected sound speed", "[eos]") {
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    REQUIRE_THAT(euler::sound_speed(right, air),
                 Catch::Matchers::WithinRel(std::sqrt(1.12), 1e-12));
}

TEST_CASE("Sod left state has the expected flux", "[eos]") {
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::IdealGas air{};
    const euler::Flux left_f = euler::physical_flux(left, air);
    REQUIRE_THAT(left_f.mass, Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(left_f.momentum, Catch::Matchers::WithinRel(1.0, 1e-12));
    REQUIRE_THAT(left_f.energy, Catch::Matchers::WithinAbs(0.0, 1e-12));
}

TEST_CASE("Sod right state has the expected flux", "[eos]") {
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    const euler::Flux right_f = euler::physical_flux(right, air);
    REQUIRE_THAT(right_f.mass, Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(right_f.momentum, Catch::Matchers::WithinRel(0.1, 1e-12));
    REQUIRE_THAT(right_f.energy, Catch::Matchers::WithinAbs(0.0, 1e-12));
}

// Sod states are at rest (u = 0), so every term containing u vanishes there.
// The moving states below exercise those terms.

TEST_CASE("Moving state converts to the expected conserved variables", "[eos]") {
    const euler::Primitive w{2.0, 3.0, 4.0};
    const euler::IdealGas air{};
    const euler::Conserved q = euler::to_conserved(w, air);
    REQUIRE_THAT(q.rho, Catch::Matchers::WithinRel(2.0, 1e-12));
    REQUIRE_THAT(q.mom, Catch::Matchers::WithinRel(6.0, 1e-12));
    REQUIRE_THAT(q.E, Catch::Matchers::WithinRel(19.0, 1e-12));
}

TEST_CASE("Moving state has the expected flux", "[eos]") {
    const euler::Primitive w{1.5, 2.0, 2.4};
    const euler::IdealGas air{};
    const euler::Flux f = euler::physical_flux(w, air);
    REQUIRE_THAT(f.mass, Catch::Matchers::WithinRel(3.0, 1e-12));
    REQUIRE_THAT(f.momentum, Catch::Matchers::WithinRel(8.4, 1e-12));
    REQUIRE_THAT(f.energy, Catch::Matchers::WithinRel(22.8, 1e-12));
}

TEST_CASE("Conserved state converts to the expected primitive variables", "[eos]") {
    const euler::Conserved q{1.5, 3.0, 9.0};
    const euler::IdealGas air{};
    const euler::Primitive w = euler::to_primitive(q, air);
    REQUIRE_THAT(w.rho, Catch::Matchers::WithinRel(1.5, 1e-12));
    REQUIRE_THAT(w.u, Catch::Matchers::WithinRel(2.0, 1e-12));
    REQUIRE_THAT(w.p, Catch::Matchers::WithinRel(2.4, 1e-12));
}

TEST_CASE("Primitive to conserved and back recovers the original state", "[eos]") {
    const euler::Primitive w{2.0, 3.0, 4.0};
    const euler::IdealGas air{};
    const euler::Primitive back = euler::to_primitive(euler::to_conserved(w, air), air);
    REQUIRE_THAT(back.rho, Catch::Matchers::WithinRel(w.rho, 1e-12));
    REQUIRE_THAT(back.u, Catch::Matchers::WithinRel(w.u, 1e-12));
    REQUIRE_THAT(back.p, Catch::Matchers::WithinRel(w.p, 1e-12));
}
