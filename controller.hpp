#pragma once

class TimestepController {
 public:
    TimestepController() = delete;

    TimestepController(double beta0, int order)
        : alpha_{1. / order - 0.75 * beta0 / order}, beta_{beta0 / order} {}

    // Determine whether previous step was successful based on error.  Calculate
    // next stepsize. If false, modify step_size. If true, modify next_step_size.
    bool success(double error, double &step_size, double &next_step_size);

 private:
    const double alpha_;
    const double beta_;
    double error_old_ = 1e-4;
    bool rejected_prev_step_ = false;
};
