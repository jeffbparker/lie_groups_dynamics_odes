#pragma once

#include <Eigen/Dense>

#include "so3.hpp"

// Represent the dynamic state of an object whose degrees of freedom live in SO(3).  E.g., the heavy
// top.
struct State {
    SO3 R;              // rotation matrix
    Eigen::Vector3d v;  // angular velocity
};

struct StateDot {
    Eigen::Vector3d v;     // angular velocity
    Eigen::Vector3d vdot;  // time derivative of velocity
};
