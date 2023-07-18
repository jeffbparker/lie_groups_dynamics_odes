#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "Eigen/Dense"
#include "adaptive_steppers.hpp"
#include "drivers.hpp"
#include "gtest/gtest.h"
#include "heavy_top.hpp"
#include "non_adaptive_steppers.hpp"
#include "so3.hpp"
#include "state.hpp"

namespace {

State initial_condition() {
    constexpr double theta0 = M_PI / 10;
    const SO3 R_initial = exp(Eigen::Vector3d{0, theta0, 0});
    const Eigen::Vector3d v_initial{0., 0., 100.};
    return {R_initial, v_initial};
}

HeavyTop create_heavy_top() {
    // Heavy top physical parameters
    constexpr double Mgl = 20;
    const Eigen::Vector3d Ib_diag{1., 1., 0.2};
    return {Mgl, Ib_diag};
}

struct LinRegressionResult {
    double intercept{};
    double slope{};
};

// Single-variable, least-squares linear regression
LinRegressionResult linear_regression(const std::vector<double>& x, const std::vector<double>& y) {
    // return on bad inputs
    if (x.size() != y.size() || x.size() == 0) return {};

    const int n = x.size();
    const double avg_x = std::accumulate(x.cbegin(), x.cend(), double{}) / n;
    const double avg_y = std::accumulate(y.cbegin(), y.cend(), double{}) / n;
    const double avg_xy = std::inner_product(x.cbegin(), x.cend(), y.cbegin(), double{}) / n;
    const double avg_xsq = std::inner_product(x.cbegin(), x.cend(), x.cbegin(), double{}) / n;

    const double slope = (avg_xy - avg_x * avg_y) / (avg_xsq - avg_x * avg_x);
    const double intercept = avg_y - slope * avg_x;
    return {intercept, slope};
}

// Assume that input data {x,y} has power-law functional form, y = C * x^p.
// Estimate p through regression.
double compute_power_exponent(const std::vector<double>& x, const std::vector<double>& y) {
    const int n = x.size();
    std::vector<double> ln_x(n);
    std::vector<double> ln_y(n);
    for (int j = 0; j < n; ++j) {
        ln_x[j] = std::log(x[j]);
        ln_y[j] = std::log(y[j]);
    }
    const LinRegressionResult result = linear_regression(ln_x, ln_y);
    return result.slope;
}

double frobenius_error(const Eigen::Matrix3d& R1, const Eigen::Matrix3d& R2) {
    return (R1 - R2).norm();
}

// Class for global shared resource
class PrecomputedSolutionEnvironment : public testing::Environment {
 public:
    void SetUp() override {
        std::cout << "Performing Global Set Up..." << std::endl;
        const State y_initial{initial_condition()};
        const HeavyTop heavy_top{create_heavy_top()};

        // ACTION
        constexpr double t_start = 0;
        constexpr double t_end = 0.4;  // 0.4

        // compute reference value with small dt.
        const double dt_ref = 6e-6;
        DriverConstantTimestep<StepperLieRK4> driver{y_initial, t_start, t_end, dt_ref, heavy_top};
        const State y_final_ref = driver.integrate();
        R_ref_ = y_final_ref.R.matrix();
        v_ref_ = y_final_ref.v.matrix();
    }

    // Reference values, computed with high precision once to be used for all tests.
    inline static Eigen::Matrix3d R_ref_ = Eigen::Matrix3d::Identity();
    inline static Eigen::Vector3d v_ref_ = Eigen::Vector3d::Zero();
};

// define a global variable
testing::Environment* const precomputed_solution_env =
    testing::AddGlobalTestEnvironment(new PrecomputedSolutionEnvironment);

}  // namespace

// Using the Heavy Top as a testbed, test the numerical properties of the Lie group integrators.

// Constant-timestep steppers: Define a fixture and Typed Test.
template <typename T>
class ConstantTimestepIntegratorTest : public ::testing::Test {};

// Associate types (constant-timestep steppers) for this Typed Test
using ConstantTimestepSteppers = testing::
    Types<StepperLieEuler, StepperLieRK2Midpoint, StepperLieRK4, StepperLieRK2CF, StepperLieRK3CF>;
TYPED_TEST_SUITE(ConstantTimestepIntegratorTest, ConstantTimestepSteppers);

// Define the tolerances for each stepper.
template <typename T>
struct ConstantTimestepToleranceTypeMap {
    static const double TOL;
};

// Default tolerance
template <typename T>
const double ConstantTimestepToleranceTypeMap<T>::TOL = 0.05;

// Tolerance for StepperLieRK3CF
template <>
const double ConstantTimestepToleranceTypeMap<StepperLieRK3CF>::TOL = 0.15;

// Define the test for each constant-timestep stepper as a heavy-top simulation.
TYPED_TEST(ConstantTimestepIntegratorTest, ErrorScalingWithDt) {
    using Stepper = TypeParam;  // convenience definition
    // SETUP
    const State y_initial{initial_condition()};
    const HeavyTop heavy_top{create_heavy_top()};

    // ACTION
    constexpr double t_start = 0;
    constexpr double t_end = 0.4;

    // loop over dt
    const std::vector<double> dt_vec{1e-5, 3e-5, 1e-4, 3e-4};
    std::vector<double> error_vec;

    for (double dt : dt_vec) {
        DriverConstantTimestep<Stepper> driver{y_initial, t_start, t_end, dt, heavy_top};
        const State y_final = driver.integrate();
        const Eigen::Matrix3d& R_final = y_final.R.matrix();
        const double error = frobenius_error(R_final, PrecomputedSolutionEnvironment::R_ref_);
        error_vec.push_back(error);
    }

    // VERIFICATION
    const double error_order = compute_power_exponent(dt_vec, error_vec);
    const double TOL = ConstantTimestepToleranceTypeMap<Stepper>::TOL;

    // asymptotically, error should scale as err ~ dt^p, with p=1 for Euler.
    EXPECT_NEAR(error_order, Stepper::EXPECTED_GLOBAL_ERROR_ORDER, TOL);
}

// Adaptive-timestep steppers: Define a fixture and Typed Test.
template <typename T>
class AdaptiveTimestepIntegratorTest : public ::testing::Test {};

// Associate types (adaptive-timestep steppers) for this Typed Test
using AdaptiveTimestepSteppers = testing::Types<StepperLieRK12, StepperLieRK23CF>;
TYPED_TEST_SUITE(AdaptiveTimestepIntegratorTest, AdaptiveTimestepSteppers);

// Define the test for each adaptive-timestep stepper as a heavy-top simulation.
TYPED_TEST(AdaptiveTimestepIntegratorTest, ErrorScalingWithTol) {
    using AdaptiveStepper = TypeParam;  // convenience definition
    // SETUP
    const State y_initial{initial_condition()};
    const HeavyTop heavy_top{create_heavy_top()};
    const double initial_dt = 1e-3;
    const double min_dt = 0;
    const double max_dt = 1e2;
    const double controller_parameter_beta0 = 0;

    // ACTION
    constexpr double t_start = 0;
    constexpr double t_end = 0.4;  // 0.4

    // loop over tol
    const std::vector<double> tol_vec{1e-6, 1e-5, 1e-4};
    std::vector<double> error_vec;

    for (double tol : tol_vec) {
        const double atol = tol;
        const double rtol = tol;
        DriverAdaptiveTimestep<AdaptiveStepper> driver{y_initial,
                                                       t_start,
                                                       t_end,
                                                       atol,
                                                       rtol,
                                                       initial_dt,
                                                       min_dt,
                                                       max_dt,
                                                       controller_parameter_beta0,
                                                       heavy_top};
        const State y_final = driver.integrate();
        const Eigen::Matrix3d& R_final = y_final.R.matrix();
        const double error = frobenius_error(R_final, PrecomputedSolutionEnvironment::R_ref_);
        error_vec.push_back(error);
    }

    // VERIFICATION
    const double error_order = compute_power_exponent(tol_vec, error_vec);
    constexpr double TOL = 0.1;

    // asymptotically, error should scale as err ~ tol^1
    constexpr int EXPECTED_ERROR_ORDER = 1;
    EXPECT_NEAR(error_order, EXPECTED_ERROR_ORDER, TOL);
}
