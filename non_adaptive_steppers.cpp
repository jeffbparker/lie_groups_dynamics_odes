#include "non_adaptive_steppers.hpp"

#include "so3.hpp"

void StepperLieEuler::step(const double t, const double dt, const RHSFunctionType& rhs_func,
                           State& y) {
    // Compute dy/dt into dydt_
    rhs_func(t, y, dydt_);

    // Update state
    y.R *= exp(dt * dydt_.v);
    y.v += dt * dydt_.vdot;
}

void StepperLieRK2Midpoint::step(const double t, const double dt, const RHSFunctionType& rhs_func,
                                 State& y) {
    // Stage 1. Compute dy/dt into dydt_
    rhs_func(t, y, dydt1_);
    y2_.R = y.R * exp(dt / 2 * y.v);
    y2_.v = y.v + dt / 2 * dydt1_.vdot;

    // Stage 2.
    rhs_func(t + dt, y2_, dydt2_);

    // Update state
    y.R *= exp(dt * dydt2_.v);
    y.v += dt * dydt2_.vdot;
}
