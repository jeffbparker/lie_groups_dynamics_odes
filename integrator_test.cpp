#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "Eigen/Dense"
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

}  // namespace

// Using the Heavy Top as a testbed, test the numerical properties of the Lie group integrators.
// The expected solution has been pre-computed separately with standard RK4 at high precision (small
// timestep).

TEST(LieGroupIntegratorTest, LieEulerErrorScalingWithDt) {
    // SETUP
    const State y_initial{initial_condition()};
    const HeavyTop heavy_top{create_heavy_top()};
    const double R22_expected = 0.948719679138258;

    // ACTION
    constexpr double t_start = 0;
    constexpr double t_end = 0.4;  // 0.4

    // loop over dt
    const std::vector<double> dt_vec{1e-6, 1e-5, 1e-4};
    std::vector<double> error_vec;

    for (double dt : dt_vec) {
        DriverConstantTimestep<StepperLieEuler> driver{y_initial, t_start, t_end, dt, heavy_top};
        const State y_final = driver.integrate();
        const double R22 = y_final.R.matrix()(2, 2);
        const double error = std::abs((R22 - R22_expected) / R22_expected);
        error_vec.push_back(error);
    }

    // VERIFICATION
    const double error_order = compute_power_exponent(dt_vec, error_vec);
    constexpr double TOL = 0.05;

    // asymptotically, error should scale as err ~ dt^p, with p=1 for Euler.
    std::cout << "error order = " << error_order << std::endl;
    constexpr int EXPECTED_ERROR_ORDER = 1;
    EXPECT_NEAR(error_order, EXPECTED_ERROR_ORDER, TOL);
}

TEST(LieGroupIntegratorTest, LieRK2ErrorScalingWithDt) {
    // SETUP
    const State y_initial{initial_condition()};
    const HeavyTop heavy_top{create_heavy_top()};
    const double R22_expected = 0.948719679138258;

    // ACTION
    constexpr double t_start = 0;
    constexpr double t_end = 0.4;  // 0.4

    // loop over dt
    const std::vector<double> dt_vec{1e-5, 1e-4, 1e-3};
    std::vector<double> error_vec;

    for (double dt : dt_vec) {
        DriverConstantTimestep<StepperLieRK2Midpoint> driver{y_initial, t_start, t_end, dt,
                                                             heavy_top};
        const State y_final = driver.integrate();
        const double R22 = y_final.R.matrix()(2, 2);
        const double error = std::abs((R22 - R22_expected) / R22_expected);
        error_vec.push_back(error);
    }

    // VERIFICATION
    const double error_order = compute_power_exponent(dt_vec, error_vec);
    constexpr double TOL = 0.05;

    // asymptotically, error should scale as err ~ dt^p, with p=2 for Lie-RK2 Midpoint.
    std::cout << "error order = " << error_order << std::endl;
    constexpr int EXPECTED_ERROR_ORDER = 2;
    EXPECT_NEAR(error_order, EXPECTED_ERROR_ORDER, TOL);
}
