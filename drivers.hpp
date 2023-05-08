#pragma once

#include <iostream>
#include <utility>

#include "aliases.hpp"
#include "state.hpp"

// Class for constant-timestep integrators
template <typename Stepper>
class DriverConstantTimestep {
 public:
    DriverConstantTimestep(State y_start, double t_start, double t_end, double desired_dt,
                           RHSFunctionType rhs_func)
        : y_{std::move(y_start)},
          t_{t_start},
          t_start_{t_start},
          t_end_{t_end},
          desired_dt_{desired_dt},
          rhs_func_{std::move(rhs_func)},
          stepper_{} {}

    State integrate() {
        const double duration = t_end_ - t_start_;
        const int num_steps = static_cast<int>(std::ceil(duration / desired_dt_));

        const double actual_dt = duration / num_steps;
        for (int j = 0; j < num_steps; ++j) {
            stepper_.step(t_, actual_dt, rhs_func_, y_);
            t_ += actual_dt;
        }
        return y_;
    }

 private:
    State y_;
    double t_;
    double t_start_;
    double t_end_;
    double desired_dt_;
    RHSFunctionType rhs_func_;
    Stepper stepper_;  // Timestepper, e.g., Euler, Runge-Kutta 4th order, etc.
};
