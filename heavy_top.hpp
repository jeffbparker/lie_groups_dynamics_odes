#pragma once

#include "Eigen/Dense"
#include "state.hpp"

// Lie Group version for function f on the RHS of the ODE dy/dt = f(t, y).
class HeavyTop {
 public:
    HeavyTop(double Mgl, const Eigen::Vector3d& Ib_diag) : Mgl_{Mgl}, Ib_diag_{Ib_diag} {}

    // physical parameters
    double Mgl_;               // mass * gravitational acceleration * length
    Eigen::Vector3d Ib_diag_;  // principal moments of inertia

    // The differential equation for the heavy top system is given by:
    //
    // dR/dt  =  R * hat(v)     R in SO(3), v in R^3, hat(v) in so(3).
    // dv/dt  =  f(R, v)    [this function is computed below in dydt.vdot]
    void operator()(double t, const State& y, StateDot& dydt);
};
