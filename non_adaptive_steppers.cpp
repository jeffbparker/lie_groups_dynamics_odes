#include "non_adaptive_steppers.hpp"

#include "so3.hpp"

void StepperLieEuler::step(const double t,
                           const double dt,
                           const RHSFunctionType& rhs_func,
                           State& y) {
    // Compute dy/dt into dydt_
    rhs_func(t, y, dydt_);

    // Update state
    y.R *= exp(dt * dydt_.v);
    y.v += dt * dydt_.vdot;
}

void StepperLieRK2Midpoint::step(const double t,
                                 const double dt,
                                 const RHSFunctionType& rhs_func,
                                 State& y) {
    // Stage 1. Compute dy/dt into dydt_
    rhs_func(t, y, dydt1_);

    // Stage 2
    y2_.R = y.R * exp(dt / 2 * y.v);
    y2_.v = y.v + dt / 2 * dydt1_.vdot;
    rhs_func(t + dt / 2, y2_, dydt2_);

    // Update state
    y.R *= exp(dt * dydt2_.v);
    y.v += dt * dydt2_.vdot;
}

void StepperLieRK4::step(const double t,
                         const double dt,
                         const RHSFunctionType& rhs_func,
                         State& y) {
    // Stage 1. Compute dy/dt into dydt_
    rhs_func(t, y, dydt1_);

    // Stage 2
    const so3_vec theta2 = dt / 2 * dydt1_.v;
    ytemp_.R = y.R * exp(theta2);
    ytemp_.v = y.v + dt / 2 * dydt1_.vdot;
    rhs_func(t + dt / 2, ytemp_, dydt2_);

    // Stage 3
    const so3_vec theta3 = dt * (dydt2_.v / 2 + dt * ad(dydt1_.v, dydt2_.v) / 8);
    ytemp_.R = y.R * exp(theta3);
    ytemp_.v = y.v + dt / 2 * dydt2_.vdot;
    rhs_func(t + dt / 2, ytemp_, dydt3_);

    // Stage 4
    const so3_vec theta4 = dt * dydt3_.v;
    ytemp_.R = y.R * exp(theta4);
    ytemp_.v = y.v + dt * dydt3_.vdot;
    rhs_func(t + dt, ytemp_, dydt4_);

    // Update state
    const so3_vec theta_end =
        dt / 6 *
        (dydt1_.v + 2 * dydt2_.v + 2 * dydt3_.v + dydt4_.v + dt * ad(dydt1_.v, dydt4_.v) / 2);
    y.R *= exp(theta_end);
    y.v += dt / 6 * (dydt1_.vdot + 2 * dydt2_.vdot + 2 * dydt3_.vdot + dydt4_.vdot);
}

void StepperLieRK2CF::step(const double t,
                           const double dt,
                           const RHSFunctionType& rhs_func,
                           State& y) {
    // Stage 1. Compute dy/dt into dydt_
    rhs_func(t, y, dydt1_);

    // Stage 2
    const so3_vec theta2 = dt / 3. * dydt1_.v;
    y2_.R = y.R * exp(theta2);
    y2_.v = y.v + dt / 3. * dydt1_.vdot;
    rhs_func(t + dt / 3., y2_, dydt2_);

    // Stage 3
    const so3_vec theta3 = 2. * dt / 3. * dydt2_.v;
    y3_.R = y.R * exp(theta3);
    y3_.v = y.v + 2. * dt / 3. * dydt2_.vdot;
    rhs_func(t + 2. * dt / 3., y3_, dydt3_);

    // Update state
    const so3_vec theta_end = dt / 2 * (dydt2_.v + dydt3_.v);
    y.R *= exp(theta_end);
    y.v += dt / 2 * (dydt2_.vdot + dydt3_.vdot);
}

void StepperLieRK3CF::step(const double t,
                           const double dt,
                           const RHSFunctionType& rhs_func,
                           State& y) {
    // Stage 1. Compute dy/dt into dydt_
    rhs_func(t, y, dydt1_);

    // Stage 2
    const so3_vec theta2 = dt / 3. * dydt1_.v;
    y2_.R = y.R * exp(theta2);
    y2_.v = y.v + dt / 3. * dydt1_.vdot;
    rhs_func(t + dt / 3., y2_, dydt2_);

    // Stage 3
    const so3_vec theta3 = 2. * dt / 3. * dydt2_.v;
    y3_.R = y.R * exp(theta3);
    y3_.v = y.v + 2. * dt / 3. * dydt2_.vdot;
    rhs_func(t + 2. * dt / 3., y3_, dydt3_);

    // Update state
    const so3_vec theta_end = dt * (-1. / 12. * dydt1_.v + 3. / 4. * dydt3_.v);
    y.R = y2_.R * exp(theta_end);
    y.v = y2_.v + dt * (-1. / 12. * dydt1_.vdot + 3. / 4. * dydt3_.vdot);
}
