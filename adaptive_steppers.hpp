#pragma once

#include <cmath>
#include <limits>
#include <utility>

#include "aliases.hpp"
#include "controller.hpp"
#include "state.hpp"

// Class to perform a single step using the RK12 method, using the embedded pair of RK2 Midpoint and
// explicit Euler.
class StepperRK12 {
 public:
    StepperRK12(double atol, double rtol, double controller_parameter_beta0)
        : atol_{atol}, rtol_{rtol}, controller_{controller_parameter_beta0, ORDER_} {}

    void step(double try_dt,
              const RHSFunctionType& rhs_func,
              const StateDot& dydt1,
              State& y,
              double& t,
              double& performed_dt,
              double& next_dt);

 private:
    double tentative_step(double t,
                          double dt,
                          const RHSFunctionType& rhs_func,
                          const State& y,
                          const StateDot& dydt,
                          State& y_out);
    static constexpr int ORDER_ = 2;
    const double atol_;
    const double rtol_;
    static constexpr int num_vel_vars_ = 3;  // number of *velocity* variables in ODE.

    State y2_;        // state at intermediate time t + h/2, h = stepsize
    StateDot dydt2_;  // time derivative at intermediate time t + h/2
    State y_out_;

    TimestepController controller_;
};
