#pragma once

struct SimConfig
{
    int nodes {20};
    int phi_iterations {500};
    double length {1.0};
    double endTime {10.0};
    double cfl {0.1};
    double velocity_bound_cond {1.0};
    double visc {0.01};
    double velocity_threshold {velocity_bound_cond * 0.000001};
    double dx {length / (nodes - 1)};
};

inline constexpr SimConfig config {};
