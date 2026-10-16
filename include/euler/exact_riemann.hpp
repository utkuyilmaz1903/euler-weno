#pragma once

#include "euler/eos.hpp"
#include "euler/state.hpp"

namespace euler {

// Velocity change across one wave of the Riemann problem when the pressure of a side goes
// from its initial value p_K to p (Toro, "Riemann Solvers and Numerical Methods for Fluid
// Dynamics", chapter 4, the function f_K):
//   p <= p_K, rarefaction:  f_K = 2 c_K / (gamma - 1) * ((p / p_K)^((gamma - 1) / (2 gamma)) - 1)
//   p >  p_K, shock:        f_K = (p - p_K) * sqrt(A_K / (p + B_K)),
//                           A_K = 2 / ((gamma + 1) rho_K),  B_K = (gamma - 1) / (gamma + 1) p_K
// The star velocity follows from either side: u* = u_L - f_L(p*) = u_R + f_R(p*).
[[nodiscard]] double pressure_function(double p, Primitive side, IdealGas gas) noexcept;

} // namespace euler
