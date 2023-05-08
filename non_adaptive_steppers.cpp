#include "non_adaptive_steppers.hpp"

void StepperEuler::step(const double t, const double dt, const RHSFunctionType& rhs_func, State& y) {
    // Compute dy/dt into dydt_
    rhs_func(t, y, dydt_);

    // Update state
    y.R *= exp(dt * dydt_.v);
    y.v += dt * dydt_.vdot;
}
