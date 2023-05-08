#include "so3.hpp"

#include <cmath>

SO3 exp(const so3_vec& x) {
    const double theta = x.norm();
    const double a = cos(theta);
    double b, c;
    constexpr double SMALL_THETA_THRESHOLD = 0.1;
    if (theta > SMALL_THETA_THRESHOLD) {  // normal Rodriguez formula
        const double inv_th = 1 / theta;
        b = sin(theta) * inv_th;
        c = (1 - a) * inv_th * inv_th;
    } else {  // series expansion for small theta to retain high accuracy
        double th_sq = theta * theta;
        b = 1 + th_sq * (-1. / 6 + th_sq * (1. / 4 + th_sq * (-1. / 5040 + th_sq * (1. / 362880))));
        c = 0.5 +
            th_sq * (-1. / 24 + th_sq * (1. / 720 + th_sq * (-1. / 40320 + th_sq * (1. / 362880))));
    }

    // for SO(3), exp(x) = a * I + b * hat(x) + c * x * x^T,
    // where hat(x) is the isomorphism from R^3 to so(3) skew symmetric matrices.
    Eigen::Matrix3d out;
    out(0, 0) = a + c * x(0) * x(0);
    out(1, 0) = b * x(2) + c * x(0) * x(1);
    out(2, 0) = -b * x(1) + c * x(0) * x(2);

    out(0, 1) = -b * x(2) + c * x(0) * x(1);
    out(1, 1) = a + c * x(1) * x(1);
    out(2, 1) = b * x(0) + c * x(1) * x(2);

    out(0, 2) = b * x(1) + c * x(0) * x(2);
    out(1, 2) = -b * x(0) + c * x(1) * x(2);
    out(2, 2) = a + c * x(2) * x(2);

    return SO3(out);
}
