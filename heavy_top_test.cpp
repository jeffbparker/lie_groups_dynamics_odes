#include "heavy_top.hpp"

#include "Eigen/Dense"
#include "gtest/gtest.h"
#include "so3.hpp"
#include "state.hpp"

using Vec3 = Eigen::Vector3d;

TEST(HeavyTopTest, TestDerivative) {
    // SETUP
    constexpr double Mgl = 20;
    const Vec3 Ib_diag{1., 1., 0.2};
    HeavyTop heavy_top{Mgl, Ib_diag};

    const SO3 R = exp(Vec3{0, 0.3, 0.4});
    const Vec3 v{0., 1., 4.};
    const State y{R, v};

    // ACTION
    const double t = 0;
    StateDot dydt;
    heavy_top(t, y, dydt);

    // VERIFICATION
    const Vec3 vdot_expected{4.37520741, 5.75310646, 0.};
    constexpr double TOL = 1e-5;
    EXPECT_TRUE(dydt.vdot.isApprox(vdot_expected, TOL));
}
