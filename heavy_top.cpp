#include "heavy_top.hpp"

#include <iostream>

#include "Eigen/Dense"

void HeavyTop::operator()(double, const State& y, StateDot& dydt) {
    using Vec3 = Eigen::Vector3d;
    dydt.v = y.v;

    // unit vector in body frame, from fixed point to center of mass of top
    static const Vec3 chi_b = Vec3::UnitZ();

    // unit vector anti-parallel to gravity, in spatial frame
    static const Vec3 gamma_s = Vec3::UnitZ();

    // compute the unit vector anti-parallel to gravity in body frame
    const Vec3 gamma_b = y.R.inverse() * gamma_s;

    // const Vec3 term1 = (Ib_diag_.cwiseProduct(y.v)).cross(y.v);
    // const Vec3 term2 = Mgl_ * gamma_b.cross(chi_b).cwiseQuotient(Ib_diag_);

    // std::cout << term2 << std::endl;
    // dydt.vdot = term1 + term2;

    dydt.vdot = ((Ib_diag_.cwiseProduct(y.v)).cross(y.v) + Mgl_ * gamma_b.cross(chi_b))
                    .cwiseQuotient(Ib_diag_);
}
