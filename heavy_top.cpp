#include "heavy_top.hpp"

#include "Eigen/Dense"

// The differential equation for the heavy top system is given by:
//
// dR/dt  =  R * hat(v)     R in SO(3), v in R^3, hat(v) in so(3).
// dv/dt  =  f(R, v)    [this function is computed below in dydt.vdot]
void HeavyTop::operator()(double, const State& y, StateDot& dydt) {
    using Vec3 = Eigen::Vector3d;
    dydt.v = y.v;

    // unit vector in body frame, from fixed point to center of mass of top
    static const Vec3 chi_b = Vec3::UnitZ();

    // unit vector anti-parallel to gravity, in spatial frame
    static const Vec3 gamma_s = Vec3::UnitZ();

    // compute the unit vector anti-parallel to gravity in body frame
    const Vec3 gamma_b = y.R.inverse() * gamma_s;

    dydt.vdot = ((Ib_diag_.cwiseProduct(y.v)).cross(y.v) + Mgl_ * gamma_b.cross(chi_b))
                    .cwiseQuotient(Ib_diag_);
}
