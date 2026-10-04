#pragma once

#include "euler/eos.hpp"
#include "euler/state.hpp"

namespace euler {

// Rusanov (local Lax-Friedrichs) flux through the face between two cells:
//   F = 1/2 (F(U_L) + F(U_R)) - 1/2 * a * (U_R - U_L),
// where a = max(|u_L| + c_L, |u_R| + c_R) is the fastest wave speed on either side.
// The second term damps zig-zags that a plain average of the two fluxes cannot.
[[nodiscard]] Flux rusanov_flux(Conserved left, Conserved right, IdealGas gas) noexcept;

} // namespace euler
