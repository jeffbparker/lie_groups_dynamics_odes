#pragma once

#include <iostream>
#include <stdexcept>
#include <utility>

#include "aliases.hpp"
#include "state.hpp"

// Class for adaptive-timestep integrators
template <typename Stepper>
class DriverAdaptiveTimestep {
 public:
    DriverAdaptiveTimestep(State y_start,
                           double t_start,
                           double t_end,
                           double atol,
                           double rtol,
                           double initial_dt,
                           double min_dt,
                           double max_dt,
                           double controller_parameter_beta0,
                           RHSFunctionType rhs_func)
        : t_{t_start},
          y_{std::move(y_start)},
          t_start_{t_start},
          t_end_{t_end},
          min_dt_{min_dt},
          max_dt_{max_dt},
          dt_{initial_dt},
          rhs_func_{std::move(rhs_func)},
          stepper_{atol, rtol, controller_parameter_beta0} {}

    ~DriverAdaptiveTimestep() {
        std::cout << "Num steps accepted: " << num_steps_accepted_
                  << "   Num steps retried: " << num_steps_retried_ << std::endl;
    }

    // Numerically integrate from t_start to t_end.
    State integrate() {
        for (int j = 0; j < max_number_of_steps_; ++j) {
            // compute derivative at beginning of step here, in case it needs to be reused.
            rhs_func_(t_, y_, dydt_);

            // Reduce step size if step would overshoot the end point
            if (t_ + 1.0001 * dt_ > t_end_) {
                dt_ = t_end_ - t_;
            }

            double performed_dt = 0;
            double next_dt = 0;
            // Perform one step of the differential equation, updating y_, t_, performed_dt, and
            // next_dt
            stepper_.step(dt_, rhs_func_, dydt_, y_, t_, performed_dt, next_dt);

            if (performed_dt == dt_) {
                ++num_steps_accepted_;
            } else {
                ++num_steps_retried_;
            }

            // check if finished.
            double EPS = 1e-10;
            if (t_ * (1 + EPS) >= t_end_) {
                return y_;
            }

            // set next timestep
            dt_ = std::min(next_dt, max_dt_);

            if (dt_ < min_dt_) {
                throw std::runtime_error(
                    "Error: dt is below the minimum allowed value. Aborting calculation...");
            }
        }
        throw std::runtime_error(
            "Error: number of steps exceeds the maximum allowed.  Aborting calculation...");
    }

 private:
    double t_;
    State y_;
    StateDot dydt_;
    const double t_start_;
    const double t_end_;
    const double min_dt_;
    const double max_dt_;
    double dt_;
    RHSFunctionType rhs_func_;
    Stepper stepper_;  // adaptive timestepper

    static constexpr int max_number_of_steps_ = 1'000'000;
    int num_steps_accepted_ = 0;
    int num_steps_retried_ = 0;
};

// Class for constant-timestep integrators
template <typename Stepper>
class DriverConstantTimestep {
 public:
    DriverConstantTimestep(
        State y_start, double t_start, double t_end, double desired_dt, RHSFunctionType rhs_func)
        : t_{t_start},
          y_{std::move(y_start)},
          t_start_{t_start},
          t_end_{t_end},
          desired_dt_{desired_dt},
          rhs_func_{std::move(rhs_func)},
          stepper_{} {}

    // Numerically integrate from t_start to t_end.
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
    double t_;
    State y_;
    const double t_start_;
    const double t_end_;
    const double desired_dt_;
    RHSFunctionType rhs_func_;
    Stepper stepper_;  // Timestepper, e.g., Euler, Runge-Kutta 4th order, etc.
};
