#include "euler/exact_riemann.hpp"

#include <cmath>

namespace euler {

double pressure_function(double p, Primitive side, IdealGas gas) noexcept {
    const double gamma = gas.gamma;

    if (p <= side.p) {
        // Rarefaction: the gas expands from side.p down to p.
        const double c = sound_speed(side, gas);
        const double exponent = (gamma - 1.0) / (2.0 * gamma);
        return 2.0 * c / (gamma - 1.0) * (std::pow(p / side.p, exponent) - 1.0);
    }

    // Shock: the gas is compressed from side.p up to p.
    const double a = 2.0 / ((gamma + 1.0) * side.rho);
    const double b = (gamma - 1.0) / (gamma + 1.0) * side.p;
    return (p - side.p) * std::sqrt(a / (p + b));
}

} // namespace euler
