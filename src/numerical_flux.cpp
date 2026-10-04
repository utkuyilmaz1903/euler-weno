#include "euler/numerical_flux.hpp"

#include "euler/eos.hpp"
#include "euler/state.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace euler {

Flux rusanov_flux(Conserved left, Conserved right, IdealGas gas) noexcept {
    const Primitive left_primitive = to_primitive(left, gas);
    const Primitive right_primitive = to_primitive(right, gas);

    const Flux f_left = physical_flux(left_primitive, gas);
    const Flux f_right = physical_flux(right_primitive, gas);

    // Fastest wave speed on each side, |u| + c, and the larger of the two.
    const double speed_left = std::abs(left_primitive.u) + sound_speed(left_primitive, gas);
    const double speed_right = std::abs(right_primitive.u) + sound_speed(right_primitive, gas);
    const double a = std::max(speed_left, speed_right);

    return {
        .mass = 0.5 * (f_left.mass + f_right.mass) - 0.5 * a * (right.rho - left.rho),
        .momentum = 0.5 * (f_left.momentum + f_right.momentum) - 0.5 * a * (right.mom - left.mom),
        .energy = 0.5 * (f_left.energy + f_right.energy) - 0.5 * a * (right.E - left.E),
    };
}

} // namespace euler
