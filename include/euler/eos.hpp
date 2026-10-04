#pragma once

#include "euler/state.hpp"

namespace euler {

// Calorically perfect (ideal) gas: p = (gamma - 1) * rho * e.
struct IdealGas {
    double gamma = 1.4;
};

[[nodiscard]] Conserved to_conserved(Primitive w, IdealGas gas) noexcept;
[[nodiscard]] Primitive to_primitive(Conserved q, IdealGas gas) noexcept;

// Speed of sound, c = sqrt(gamma * p / rho).
[[nodiscard]] double sound_speed(Primitive w, IdealGas gas) noexcept;

// Physical flux of a single state, F(U) = (rho*u, rho*u^2 + p, u*(E + p)).
// For the flux through a face between two cells, see numerical_flux.hpp.
[[nodiscard]] Flux physical_flux(Primitive w, IdealGas gas) noexcept;

} // namespace euler
