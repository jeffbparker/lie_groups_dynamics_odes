#pragma once

class TimestepController {
 public:
    TimestepController() = delete;

    TimestepController(double beta0, int order)
        : alpha_{1. / order - 0.75 * beta0 / order}, beta_{beta0 / order} {}

    // Determine whether previous step was successful based on error.  Calculate
    // next timestep. If false, modify dt. If true, modify next_dt.
    bool success(double error, double& dt, double& next_dt);

 private:
    const double alpha_;
    const double beta_;
    double error_old_ = 1e-4;
    bool rejected_prev_step_ = false;
};
