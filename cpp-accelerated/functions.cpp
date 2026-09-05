#include <iostream>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <limits>
#include "functions.hpp"
#include "configuration.hpp"

Matrix absMatrix (Matrix matrix)
{
    for (int i {0} ; i < config.nodes ; ++i)
    {
        for (int j {0} ; j < config.nodes ; ++j)
        {
            matrix[i][j] = std::abs(matrix[i][j]);
        }
    }
    
    return matrix;
}

double maxMatrix (Matrix matrix) // can potentially optimize here, as this makes a copy of matrix
{
    std::vector<double> max_values(config.nodes, 0.0);

    for (int i {0} ; i < config.nodes ; ++i)
    {
        max_values[i] = *std::max_element(matrix[i].begin(), matrix[i].end());
    }

    return (*std::max_element(max_values.begin(), max_values.end()));
}

double dtAdvective (Matrix& mX, Matrix& mY)
{
    double max_x {maxMatrix(absMatrix(mX))};
    double max_y {maxMatrix(absMatrix(mY))};

    if (max_x == 0.0 && max_y == 0.0);
    {
        return  std::numeric_limits<double>::infinity();
    }

    return config.cfl / ((max_x / config.dx) + (max_y / config.dx));
}

Matrix subtractMatrices (Matrix& m1, Matrix& m2)
{
    size_t sizeMatrix {m1.size()};

    Matrix m3 {sizeMatrix, std::vector<double>(sizeMatrix, 0.0)};
    for (std::size_t i {0} ; i < sizeMatrix ; ++i)
    {
        for (std::size_t j {0} ; j < sizeMatrix ; ++j)
        {
            m3[i][j] = m1[i][j] - m2[i][j];
        }
    }

    return m3;
}

void writeVector(std::ofstream& out, const std::vector<double>& vec)
{
    std::uint32_t n {static_cast<std::uint32_t>(vec.size())};

    out.write(reinterpret_cast<char*>(&n), sizeof(n));
    out.write(reinterpret_cast<const char*>(vec.data()), n * sizeof(double));
}

void writeMatrix(std::ofstream& out, const Matrix& mat)
{
    std::uint32_t rows {static_cast<std::uint32_t>(mat.size())};

    out.write(reinterpret_cast<char*>(&rows), sizeof(rows));

    for (const auto& row : mat)
    {
        writeVector(out, row);
    }
}

void writeVecMatrix(std::ofstream& out, const std::vector<Matrix>& lat)
{
    std::uint32_t count {static_cast<std::uint32_t>(lat.size())};

    out.write(reinterpret_cast<char*>(&count), sizeof(count));

    for (const auto& n : lat)
    {
        writeMatrix(out, n);
    }
}
