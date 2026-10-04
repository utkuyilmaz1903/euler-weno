#include "euler/eos.hpp"

#include <cmath>

namespace euler {

Conserved to_conserved(Primitive w, IdealGas gas) noexcept {
    const double kinetic = 0.5 * w.rho * w.u * w.u;
    return {
        .rho = w.rho,
        .mom = w.rho * w.u,
        .E = w.p / (gas.gamma - 1.0) + kinetic,
    };
}

Primitive to_primitive(Conserved q, IdealGas gas) noexcept {
    const double u = q.mom / q.rho;
    const double kinetic = 0.5 * q.rho * u * u;
    return {
        .rho = q.rho,
        .u = u,
        .p = (gas.gamma - 1.0) * (q.E - kinetic),
    };
}

double sound_speed(Primitive w, IdealGas gas) noexcept {
    return std::sqrt(gas.gamma * w.p / w.rho);
}

Flux physical_flux(Primitive w, IdealGas gas) noexcept {
    const double E = to_conserved(w, gas).E;
    return {
        .mass = w.rho * w.u,
        .momentum = w.rho * w.u * w.u + w.p,
        .energy = w.u * (E + w.p),
    };
}

} // namespace euler
