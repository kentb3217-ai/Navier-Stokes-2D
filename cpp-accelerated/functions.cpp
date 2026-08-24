#include <iostream>
#include <algorithm>
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

    return (config.cfl * (std::min((config.dx / max_x), (config.dx / max_y))));
}
