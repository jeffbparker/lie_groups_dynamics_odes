#include "so3.hpp"

#include <iostream>

#include "Eigen/Dense"
#include "gtest/gtest.h"
#include "unsupported/Eigen/MatrixFunctions"

using Vec3 = Eigen::Vector3d;

namespace {

Eigen::Matrix3d eigen_exp(const Vec3& v) {
    Eigen::Matrix3d v_hat;
    v_hat << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
    return v_hat.exp();
}

}  // namespace

TEST(SO3Test, TestExpSmallV) {
    // test for small |v|
    // SETUP
    const so3_vec v{1e-4, 2.1e-5, 1.7e-5};

    // ACTION
    const SO3 exp_v = exp(v);

    // VERIFICATION
    Eigen::Matrix3d exp_v_expected = eigen_exp(v);

    constexpr double TOL = 1e-5;
    EXPECT_TRUE(exp_v.matrix().isApprox(exp_v_expected, TOL));
}

TEST(SO3Test, TestExpNormalV) {
    // test for non-small |v|
    const so3_vec v{2.3, 1.9, 0.3};

    // ACTION
    const SO3 exp_v = exp(v);

    // VERIFICATION
    Eigen::Matrix3d exp_v_expected = eigen_exp(v);
    std::cout << exp_v.matrix() << std::endl;

    constexpr double TOL = 1e-5;
    EXPECT_TRUE(exp_v.matrix().isApprox(exp_v_expected, TOL));
}
