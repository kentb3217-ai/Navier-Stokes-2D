#pragma once
#include <iostream>
#include <vector>
#include <cstdint>
#include <fstream>
#include "configuration.hpp"

using Matrix = std::vector<std::vector<double>>;

Matrix absMatrix (Matrix matrix);
double maxMatrix (Matrix matrix);
double dtAdvective (Matrix& m1, Matrix& m2);
Matrix subtractMatrices (Matrix& m1, Matrix& m2);
void writeVector(std::ofstream& out, const std::vector<double>& vec);
void writeMatrix(std::ofstream& out, const Matrix& mat);
void writeVecMatrix(std::ofstream& out, const std::vector<Matrix>& lat);
