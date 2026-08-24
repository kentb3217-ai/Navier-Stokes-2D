#pragma once
#include <iostream>
#include <vector>
#include "configuration.hpp"

using Matrix = std::vector<std::vector<double>>;

Matrix absMatrix (Matrix matrix);
double maxMatrix (Matrix matrix);
double dtAdvective (Matrix m1, Matrix m2);
