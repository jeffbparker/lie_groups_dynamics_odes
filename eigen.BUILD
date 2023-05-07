# Description:
#   Eigen is a C++ template library for linear algebra: vectors,
#   matrices, and related algorithms.

EIGEN_FILES = [
    "Eigen/**",
]

cc_library(
    name = "eigen",
    hdrs = glob(EIGEN_FILES),
    includes = ["."],
    visibility = ["//visibility:public"],
)
