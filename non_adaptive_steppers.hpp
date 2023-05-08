#pragma once

#include <cmath>

#include "Eigen/Dense"
#include "aliases.hpp"
#include "state.hpp"

// 1st-order explicit, Lie Euler
class StepperLieEuler {
 public:
    StepperLieEuler() = default;

    /**
     * @brief Given state y at time t, perform step and update y to time t + dt
     *
     * @param t - Time at beginning of step
     * @param dt -timestep
     * @param rhs_func - Function to compute time derivatives, dy/dt = f(t, y)
     * @param y - On entry, state at time t.  On exit, state at time t + dt
     */
    void step(double t, double dt, const RHSFunctionType& rhs_func, State& y);

 private:
    StateDot dydt_;  // time derivatives at t
};

// 2nd-order, Lie-RK2-midpoint
class StepperLieRK2Midpoint {
 public:
    StepperLieRK2Midpoint() = default;

    void step(double t, double dt, const RHSFunctionType& rhs_func, State& y);

 private:
    State y2_;        // intermediate state at t + dt/2
    StateDot dydt1_;  // time derivative at t
    StateDot dydt2_;  // time derivative at t + dt/2
};
