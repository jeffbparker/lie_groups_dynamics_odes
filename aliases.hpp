#pragma once

#include <Eigen/Dense>
#include <functional>

#include "state.hpp"

using Vec3 = Eigen::Vector3d;

// type for function f on the RHS of the ODE dy/dt = f(t, y).
using RHSFunctionType = std::function<void(double, const State&, StateDot&)>;
