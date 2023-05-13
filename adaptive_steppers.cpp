#include "adaptive_steppers.hpp"

#include <cmath>
#include <stdexcept>

void AdaptiveStepper::step(const double try_dt,
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

double StepperLieRK12::tentative_step(const double t,
                                      const double dt,
                                      const RHSFunctionType& rhs_func,
                                      const State& y,
                                      const StateDot& dydt1,
                                      State& y_out) {
    // Stage 1 precomputed into dydt1

    // Stage 2
    y2_.R = y.R * exp(dt / 2 * dydt1.v);
    y2_.v = y.v + dt / 2 * dydt1.vdot;
    rhs_func(t + dt / 2, y2_, dydt2_);

    // Compute y_out
    y_out.R = y.R * exp(dt * dydt2_.v);
    y_out.v = y.v + dt * dydt2_.vdot;

    // Compute the velocity error term, difference in 2rd-order & 1st-order results
    const decltype(dydt1.vdot) vel_error = dt * (-dydt1.vdot + dydt2_.vdot);

    // Compute the scalar error. Iterate over velocity variables in y.v to accumulate error
    double scalar_error = 0;
    for (int j = 0; j < num_vel_vars_; ++j) {
        const double error_scale = atol_ + rtol_ * std::max(std::abs(y.v[j]), std::abs(y_out.v[j]));
        scalar_error += vel_error[j] * vel_error[j] / (error_scale * error_scale);
    }
    return scalar_error;
}

double StepperLieRK23CF::tentative_step(const double t,
                                        const double dt,
                                        const RHSFunctionType& rhs_func,
                                        const State& y,
                                        const StateDot& dydt1,
                                        State& y_out) {
    // Stage 1 precomputed into dydt1

    // Stage 2
    const so3_vec theta2 = dt / 3. * dydt1.v;
    y2_.R = y.R * exp(theta2);
    y2_.v = y.v + dt / 3. * dydt1.vdot;
    rhs_func(t + dt / 3., y2_, dydt2_);

    // Stage 3
    const so3_vec theta3 = 2. * dt / 3. * dydt2_.v;
    y3_.R = y.R * exp(theta3);
    y3_.v = y.v + 2. * dt / 3. * dydt2_.vdot;
    rhs_func(t + 2. * dt / 3., y3_, dydt3_);

    // Compute y_out
    const so3_vec theta_end = dt * (-1. / 12. * dydt1.v + 3. / 4. * dydt3_.v);
    y_out.R = y2_.R * exp(theta_end);
    y_out.v = y2_.v + dt * (-1. / 12. * dydt1.vdot + 3. / 4. * dydt3_.vdot);

    // State y_2ndorder;
    // y_2ndorder.v = y.v + dt / 2 * (dydt2_.vdot + dydt3_.vdot);
    // const decltype(dydt1.vdot) vel_error = y_out.v - y_2ndorder.v;

    // Compute the velocity error term, difference in 3rd-order & 2nd-order results
    const decltype(dydt1.vdot) vel_error = dt * (dydt1.vdot / 4 - dydt2_.vdot / 2 + dydt3_.vdot / 4);

    // Compute the scalar error. Iterate over velocity variables in y.v to accumulate error
    double scalar_error = 0;
    for (int j = 0; j < num_vel_vars_; ++j) {
        const double error_scale = atol_ + rtol_ * std::max(std::abs(y.v[j]), std::abs(y_out.v[j]));
        scalar_error += vel_error[j] * vel_error[j] / (error_scale * error_scale);
    }
    return scalar_error;
}
