// Work in progress

#include <iostream>
#include <vector>
#include "configuration.hpp"
#include "functions.hpp"
#include "functions.cpp"

int main()
{
    double counter {0.0};
    double nodes {static_cast<double>(config.nodes)}; // to prevent narrowing conversions

    // velocity of the fluid in x and y directions respectively
    Matrix u_x {config.nodes, std::vector<double>{nodes, 0.0}};
    Matrix u_y {config.nodes, std::vector<double>{nodes, 0.0}};

    // Boundary conditions
    std::vector<double> botx {std::vector<double>{nodes, 0.0}};
    std::vector<double> boty {std::vector<double>{nodes, 0.0}};
    std::vector<double> upx {std::vector<double>{nodes, 0.0}};
    std::vector<double> upy {std::vector<double>{nodes, 0.0}};
    std::vector<double> leftx {std::vector<double>{nodes, 0.0}};
    std::vector<double> lefty {std::vector<double>{nodes, 0.0}};
    std::vector<double> rightx {std::vector<double>{nodes, 0.0}};
    std::vector<double> righty {std::vector<double>{nodes, 0.0}};

    u_x[0] = botx;
    u_y[0] = boty;
    u_x[config.nodes - 1] = upx;
    u_y[config.nodes - 1] = upy;

    for (int i {0} ; i < config.nodes ; ++i)
    {
        u_x[i][0] = leftx[i];
        u_y[i][0] = lefty[i];
    }

    for (int i {0} ; i < config.nodes ; ++i)
    {
        u_x[i][config.nodes - 1] = rightx[i];
        u_y[i][config.nodes - 1] = righty[i];
    }

    // Initial delta t calculation
    double advective {dtAdvective(u_x, u_y)};
    double diffusive {(config.dx*config.dx) / (4 * config.visc)};
    double dt {std::min(advective, diffusive)};

    // WHILE LOOP
    double counter {0.0};

    // Initializing partial derivatives, laplacians, auxiliary fields, phi, gradient phi, and divergence of the auxiliary field
    Matrix phi {config.nodes, std::vector<double>{nodes, 0.0}};

    Matrix aux_X {config.nodes, std::vector<double>{nodes, 0.0}};
    Matrix aux_Y {config.nodes, std::vector<double>{nodes, 0.0}};

    Matrix gradPhi_X {config.nodes, std::vector<double>{nodes, 0.0}};
    Matrix gradPhi_Y {config.nodes, std::vector<double>{nodes, 0.0}};

    Matrix divAuxField {config.nodes, std::vector<double>{nodes, 0.0}};

    Matrix dux_dx {config.nodes, std::vector<double>{nodes, 0.0}};
    Matrix dux_dy {config.nodes, std::vector<double>{nodes, 0.0}};

    Matrix duy_dx {config.nodes, std::vector<double>{nodes, 0.0}};
    Matrix duy_dy {config.nodes, std::vector<double>{nodes, 0.0}};

    Matrix dd_ux {config.nodes, std::vector<double>{nodes, 0.0}};
    Matrix dd_uy {config.nodes, std::vector<double>{nodes, 0.0}};

    while (counter < config.endTime)
    {
        for (int i {1} ; i < (config.nodes - 1) ; ++i)
        {
            for (int j {1} ; j < (config.nodes - 1) ; ++j)
            {
                
            }
        }
    }

    return 0;
}
