#pragma once

namespace euler {

// Primitive variables: what you would measure in the flow.
struct Primitive {
    double rho; // density
    double u;   // velocity
    double p;   // pressure
};

// Conserved variables: what the finite-volume scheme evolves.
struct Conserved {
    double rho; // density
    double mom; // momentum, rho * u
    double E;   // total energy per unit volume
};

// Flux through a face: what crosses it per unit time and unit area.
// A separate type from Conserved so the two cannot be mixed up by accident.
struct Flux {
    double mass;     // rho * u
    double momentum; // rho * u^2 + p
    double energy;   // u * (E + p)
};

} // namespace euler
