load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "state",
    hdrs = ["state.hpp"],
    deps = [
        ":so3",
        "@eigen",
    ],
)

cc_library(
    name = "aliases",
    hdrs = ["aliases.hpp"],
    deps = [
        ":state",
        "@eigen",
    ],
)

cc_library(
    name = "so3",
    srcs = ["so3.cpp"],
    hdrs = ["so3.hpp"],
    deps = [
        "@eigen",
    ],
)

cc_test(
    name = "so3_test",
    srcs = ["so3_test.cpp"],
    deps = [
        ":so3",
        "@gtest//:gtest_main",
    ],
)

cc_library(
    name = "controller",
    srcs = ["controller.cpp"],
    hdrs = ["controller.hpp"],
)

cc_library(
    name = "non_adaptive_steppers",
    srcs = ["non_adaptive_steppers.cpp"],
    hdrs = ["non_adaptive_steppers.hpp"],
    deps = [
        ":aliases",
        ":so3",
        ":state",
        "@eigen",
    ],
)

cc_library(
    name = "drivers",
    hdrs = ["drivers.hpp"],
    deps = [
        ":aliases",
        ":state",
    ],
)

cc_library(
    name = "heavy_top",
    srcs = ["heavy_top.cpp"],
    hdrs = ["heavy_top.hpp"],
    deps = [
        ":state",
    ],
)

cc_test(
    name = "heavy_top_test",
    srcs = ["heavy_top_test.cpp"],
    deps = [
        ":heavy_top",
        ":so3",
        ":state",
        "@eigen",
        "@gtest//:gtest_main",
    ],
)

cc_test(
    name = "integrator_test",
    size = "small",
    srcs = ["integrator_test.cpp"],
    deps = [
        ":drivers",
        ":heavy_top",
        ":non_adaptive_steppers",
        ":so3",
        ":state",
        "@gtest//:gtest_main",
    ],
)
