#pragma once

#include <Eigen/Dense>

// Basic operations for Lie group SO(3) needed for demonstration purposes

// Class for Lie group elements SO(3)
class SO3 {
 public:
    // Constructor, assuming R is a valid SO(3); i.e. does not check it is a valid rotation matrix.
    explicit SO3(const Eigen::Matrix3d& mat) : mat_{std::move(mat)} {}

    // Group product of two elements (matrix multiplication)
    SO3 operator*(const SO3& other) const { return SO3(mat_ * other.mat_); }

    // Update group element due to right-product with another element
    SO3& operator*=(const SO3& other) {
        mat_ *= other.mat_;
        return *this;
    }

    // Product of group element and vector
    Eigen::Vector3d operator*(const Eigen::Vector3d& v) const { return mat_ * v; }

    // Inverse element
    SO3 inverse() const { return SO3{mat_.transpose()}; }

    // Get the underlying matrix
    const Eigen::Matrix3d& matrix() const { return mat_; }
    Eigen::Matrix3d& matrix() { return mat_; }

 private:
    Eigen::Matrix3d mat_;
};

// alias for the R^3 vector representation isomorphic to Lie algebra so(3).
// All computations are done with the R^3 representation rather than with skew-symmetric matrices.
using so3_vec = Eigen::Vector3d;

// exp: so(3) -> SO(3)
SO3 exp(const so3_vec& x);

// NOT IMPLEMENTED YET
// // adjoint representation.  ad(x, y) = [x, y] = xy - yx
// so3_vec ad(const so3_vec& x, const so3_vec& y);

// // inverse of dexp.
// so3_vec dexpinv(const so3_vec& x, const so3_vec& y);
