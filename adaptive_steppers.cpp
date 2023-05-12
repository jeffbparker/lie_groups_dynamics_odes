#include "adaptive_steppers.hpp"

#include <cmath>
#include <stdexcept>

void StepperRK12::step(const double try_dt,
                       const RHSFunctionType& rhs_func,
                       const StateDot& dydt1,
                       State& y,
                       double& t,
                       double& performed_dt,
                       double& next_dt) {
    double dt = try_dt;  // initialize dt to incoming value
    while (true) {
        double scalar_error = tentative_step(t, dt, rhs_func, y, dydt1, y_out_);
        if (controller_.success(scalar_error, dt, next_dt)) {
            break;
        }

        if (std::abs(dt) <= std::abs(t) * std::numeric_limits<double>::epsilon()) {
            throw std::runtime_error("Error: timestep underflow.  Aborting calculation...");
        }
    }

    std::swap(y, y_out_);
    performed_dt = dt;
    t += dt;
}

double StepperRK12::tentative_step(const double t,
                                   const double dt,
                                   const RHSFunctionType& rhs_func,
                                   const State& y,
                                   const StateDot& dydt,
                                   State& y_out) {
    // stage 1 precomputed into dydt
    // stage 2
    y2_.R = y.R * exp(dt / 2 * dydt.v);
    y2_.v = y.v + dt / 2 * dydt.vdot;
    rhs_func(t + dt / 2, y2_, dydt2_);

    // Compute y_out
    y_out.R = y.R * exp(dt * dydt2_.v);
    y_out.v = y.v + dt * dydt2_.vdot;

    // Compute the scalar error. Iterate over velocity variables in y.v to accumulate error
    const decltype(y.v) vel_error = dt * (-dydt.vdot + dydt2_.vdot);
    double scalar_error = 0;
    for (int j = 0; j < num_vel_vars_; ++j) {
        const double error_scale = atol_ + rtol_ * std::max(std::abs(y.v[j]), std::abs(y_out.v[j]));
        scalar_error += vel_error[j] * vel_error[j] / (error_scale * error_scale);
    }
    return scalar_error;
}
