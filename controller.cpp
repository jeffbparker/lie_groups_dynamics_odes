#include "controller.hpp"

#include <algorithm>
#include <cmath>

bool TimestepController::success(const double error, double &step_size, double &next_step_size) {
    constexpr double SAFETY_FACTOR = 0.9;

    // Min and max scale factors by which the step size can change in one step.
    constexpr double MIN_SCALE_FACTOR = 0.2;
    constexpr double MAX_SCALE_FACTOR = 5.0;

    double scale_factor = 0.0;
    if (error <= 1.0) {  // success; compute next_step_size
        if (error == 0.0) {
            scale_factor = MAX_SCALE_FACTOR;
        } else {
            scale_factor = SAFETY_FACTOR * std::pow(error, -alpha_) * std::pow(error_old_, beta_);
            scale_factor = std::clamp(scale_factor, MIN_SCALE_FACTOR, MAX_SCALE_FACTOR);
        }

        // If the prior step attempt was rejected, don't let the step size increase
        if (rejected_prev_step_) {
            next_step_size = step_size * std::min(scale_factor, 1.0);
        } else {
            next_step_size = step_size * scale_factor;
        }

        // store the error for the next call
        constexpr double MIN_STORED_ERROR = 1.0e-4;
        error_old_ = std::max(error, MIN_STORED_ERROR);
        rejected_prev_step_ = false;
        return true;
    } else {  // step failed.  Compute the reduced stepsize for next step attempt.
        scale_factor = SAFETY_FACTOR * std::pow(error, -alpha_);
        scale_factor = std::max(scale_factor, MIN_SCALE_FACTOR);
        step_size *= scale_factor;
        rejected_prev_step_ = true;
        return false;
    }
}
