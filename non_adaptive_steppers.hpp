#pragma once

#include <cmath>

#include "Eigen/Dense"
#include "aliases.hpp"
#include "state.hpp"

// 1st-order explicit, Lie Euler
class StepperLieEuler {
 public:
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
    State y2_;        // intermediate state at t + dt/2
    StateDot dydt1_;  // time derivative at t
    StateDot dydt2_;  // time derivative at t + dt/2
};

// 4th-order, Lie-RK4, graded free algebra minimizing the number of commutators
class StepperLieRK4 {
 public:
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
    // intermediate-stage states and derivative storage
    State ytemp_;

    StateDot dydt1_;
    StateDot dydt2_;
    StateDot dydt3_;
    StateDot dydt4_;
};

// 2nd-order Lie-RK2CF (Commutator free), part of 3-2 adaptive method
class StepperLieRK2CF {
 public:
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
    // intermediate-stage states and derivative storage
    State y2_;
    State y3_;

    StateDot dydt1_;
    StateDot dydt2_;
    StateDot dydt3_;
};

// 3rd-order Lie-RK3CF (Commutator free), part of 3-2 adaptive method
class StepperLieRK3CF {
 public:
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
    // intermediate-stage states and derivative storage
    State y2_;
    State y3_;

    StateDot dydt1_;
    StateDot dydt2_;
    StateDot dydt3_;
};
