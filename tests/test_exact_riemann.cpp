#include "euler/eos.hpp"
#include "euler/exact_riemann.hpp"
#include "euler/state.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Sod: left (1, 0, 1), right (0.125, 0, 0.1), air.

TEST_CASE("No pressure change means no velocity change on either side", "[exact]") {
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    const double f_left = euler::pressure_function(1.0, left, air);   // p = p_L
    const double f_right = euler::pressure_function(0.1, right, air); // p = p_R
    REQUIRE_THAT(f_left, Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(f_right, Catch::Matchers::WithinAbs(0.0, 1e-12));
}

// p = 0.2 is below p_L (left gas expands: rarefaction branch) and above p_R (right gas is
// compressed: shock branch), so one call exercises each branch.
TEST_CASE("Pressure function matches known values on both branches", "[exact]") {
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    const double f_left = euler::pressure_function(0.2, left, air);
    const double f_right = euler::pressure_function(0.2, right, air);
    REQUIRE_THAT(f_left, Catch::Matchers::WithinAbs(-1.215178, 1e-6));
    REQUIRE_THAT(f_right, Catch::Matchers::WithinAbs(0.554700, 1e-6));
}

// At p* both sides give the same star velocity, u_L - f_L = u_R + f_R. With u_L = u_R = 0
// that means f_L + f_R = 0; each one alone is about -/+0.9275, not zero.
TEST_CASE("At the Sod star pressure the two sides balance", "[exact]") {
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    const double p_star = 0.30313; // Toro, table 4.2
    const double f_left = euler::pressure_function(p_star, left, air);
    const double f_right = euler::pressure_function(p_star, right, air);
    REQUIRE_THAT(f_left + f_right, Catch::Matchers::WithinAbs(0.0, 1e-5));
}

// The derivative formulas were taken from Toro without derivation, and Newton would still
// converge with a slightly wrong slope (only more slowly), so check them directly against the
// slope measured from pressure_function, (f(p + h) - f(p)) / h. With h = 1e-6 the two agree to
// about 2e-6 relative: a larger h measures over a more curved stretch, a much smaller one loses
// digits when subtracting two nearly equal values.
TEST_CASE("Pressure function derivative matches the measured slope on both branches", "[exact]") {
    const euler::Primitive left{1.0, 0.0, 1.0};
    const euler::Primitive right{0.125, 0.0, 0.1};
    const euler::IdealGas air{};
    const double p = 0.2; // left gas on the rarefaction branch, right gas on the shock branch
    const double h = 1e-6;

    const double f_left = euler::pressure_function(p, left, air);
    const double f_left_h = euler::pressure_function(p + h, left, air);
    const double f_right = euler::pressure_function(p, right, air);
    const double f_right_h = euler::pressure_function(p + h, right, air);
    const double measured_slope_left = (f_left_h - f_left) / h;
    const double measured_slope_right = (f_right_h - f_right) / h;

    const double derivative_left = euler::pressure_function_derivative(p, left, air);
    const double derivative_right = euler::pressure_function_derivative(p, right, air);
    REQUIRE_THAT(derivative_left, Catch::Matchers::WithinRel(measured_slope_left, 1e-5));
    REQUIRE_THAT(derivative_right, Catch::Matchers::WithinRel(measured_slope_right, 1e-5));
}