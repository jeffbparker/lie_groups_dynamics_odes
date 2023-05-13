#pragma once

#include <cmath>
#include <limits>
#include <utility>

#include "aliases.hpp"
#include "controller.hpp"
#include "state.hpp"

// Abstract class for performing steps using an adaptive-timestep method.
class AdaptiveStepper {
 public:
    /**
     * @brief Construct a new Adaptive Stepper object
     *
     * @param atol - absolute tolerance
     * @param rtol - relative tolerance
     * @param controller_parameter_beta0 - determines the level of PI in the timestep controller
     * @param order  - the *local* error of the *lower-order* method in the embedded pair. E.g., 5
     * for RK45
     */
    AdaptiveStepper(double atol, double rtol, double controller_parameter_beta0, int order)
        : atol_{atol}, rtol_{rtol}, controller_{controller_parameter_beta0, order} {}
    virtual ~AdaptiveStepper() = default;

    /**
     * @brief - Perform one step.  Will retry steps as needed to ensure error is small enough.
     *
     * @param[in] try_dt - Timestep to be tried first
     * @param[in] rhs_func - the RHS of the ODE, i.e. the time derivative
     * @param[in] dydt1 - The time derivative at time t computed from y, which has been precomputed
     * @param[inout] y - On entry, state at time t.  On exit, state at time t + performed_dt
     * @param[inout] t - On entry, time at beginning of step.  On exit, time at end of step
     * @param[out] performed_dt - The timestep of the accepted step.
     * @param[out] next_dt - The timestep computed by the timestep controller for the next step
     */
    void step(double try_dt,
              const RHSFunctionType& rhs_func,
              const StateDot& dydt1,
              State& y,
              double& t,
              double& performed_dt,
              double& next_dt);

    const double atol_;
    const double rtol_;
    static constexpr int num_vel_vars_ = 3;  // number of *velocity* variables in ODE.
    State y_out_;
    TimestepController controller_;

 private:
    /**
     * @brief  - Perform one tentative step and return the error estimate.  Does not retry steps.
     *
     * @param[in] t - time at beginning of step
     * @param[in] dt - Timestep
     * @param[in] rhs_func - the RHS of the ODE, i.e. the time derivative
     * @param[in] y - State at time t, prior to step
     * @param[in] dydt1 - The time derivative at time t computed from y, which has been precomputed
     * @param[out] y_out - State after the step at time t + dt, using higher-order method
     * @return double - error estimate
     */
    virtual double tentative_step(double t,
                                  double dt,
                                  const RHSFunctionType& rhs_func,
                                  const State& y,
                                  const StateDot& dydt1,
                                  State& y_out) = 0;
};

// Class to perform a single step using the RK12 method, using the embedded pair of RK2 Midpoint and
// explicit Euler.
class StepperLieRK12 final : public AdaptiveStepper {
 public:
    /**
     * @brief Construct a new StepperLieRK12 object
     *
     * @param[in] atol - absolute tolerance
     * @param[in] rtol - relative tolerance
     * @param[in] controller_parameter_beta0 - determines the level of PI in the timestep controller
     */
    StepperLieRK12(double atol, double rtol, double controller_parameter_beta0)
        : AdaptiveStepper{atol, rtol, controller_parameter_beta0, ORDER_} {}
    ~StepperLieRK12() {}

 private:
    /**
     * @brief  - Perform one tentative step and return the error estimate from the RK12 embedded
     * pair.  Does not retry steps.
     *
     * @param[in] t - time at beginning of step
     * @param[in] dt - Timestep
     * @param[in] rhs_func - the RHS of the ODE, i.e. the time derivative
     * @param[in] y - State at time t, prior to step
     * @param[in] dydt1 - The time derivative at time t computed from y, which has been precomputed
     * @param[out] y_out - State after the step at time t + dt, using higher-order method
     * @return double - error estimate
     */
    double tentative_step(double t,
                          double dt,
                          const RHSFunctionType& rhs_func,
                          const State& y,
                          const StateDot& dydt1,
                          State& y_out) override;
    static constexpr int ORDER_ = 2;

    State y2_;        // state at intermediate time t + h/2, h = stepsize
    StateDot dydt2_;  // time derivative at intermediate time t + h/2
};

// Class to perform a single step using the RK23 Commutator Free method
class StepperLieRK23CF final : public AdaptiveStepper {
 public:
    /**
     * @brief Construct a new StepperLieRK23 object
     *
     * @param[in] atol - absolute tolerance
     * @param[in] rtol - relative tolerance
     * @param[in] controller_parameter_beta0 - determines the level of PI in the timestep controller
     */
    StepperLieRK23CF(double atol, double rtol, double controller_parameter_beta0)
        : AdaptiveStepper{atol, rtol, controller_parameter_beta0, ORDER_} {}

 private:
    /**
     * @brief  - Perform one tentative step and return the error estimate from the RK23 embedded
     * pair.  Does not retry steps.
     *
     * @param[in] t - time at beginning of step
     * @param[in] dt - Timestep
     * @param[in] rhs_func - the RHS of the ODE, i.e. the time derivative
     * @param[in] y - State at time t, prior to step
     * @param[in] dydt1 - The time derivative at time t computed from y, which has been precomputed
     * @param[out] y_out - State after the step at time t + dt, using higher-order method
     * @return double - error estimate
     */
    double tentative_step(double t,
                          double dt,
                          const RHSFunctionType& rhs_func,
                          const State& y,
                          const StateDot& dydt1,
                          State& y_out) override;
    static constexpr int ORDER_ = 3;

    // intermediate-stage states and derivative storage
    State y2_;
    State y3_;

    StateDot dydt2_;
    StateDot dydt3_;
};
