// Work in progress
#include <limits>
#include <iostream>
#include <vector>
#include "configuration.hpp"
#include "functions.hpp"
#include "functions.cpp"

// Note: Go through and check if its matrix[nodes] or matrix[nodes - 1] for last col/row

int main()
{
    double nodes {static_cast<double>(config.nodes)};

    // velocity of the fluid in x and y directions respectively
    Matrix u_x {config.nodes, std::vector<double>(config.nodes, 0.0)};
    Matrix u_y {config.nodes, std::vector<double>(config.nodes, 0.0)};

    // Boundary conditions
    std::vector<double> botx {std::vector<double>(config.nodes, 0.0)};
    std::vector<double> boty {std::vector<double>(config.nodes, 0.0)};
    std::vector<double> upx {std::vector<double>(config.nodes, 0.0)};
    std::vector<double> upy {std::vector<double>(config.nodes, 0.0)};
    std::vector<double> leftx {std::vector<double>(config.nodes, 0.0)};
    std::vector<double> lefty {std::vector<double>(config.nodes, 0.0)};
    std::vector<double> rightx {std::vector<double>(config.nodes, 0.0)};
    std::vector<double> righty {std::vector<double>(config.nodes, 0.0)};

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
    Matrix phi {config.nodes, std::vector<double>(config.nodes, 0.0)};

    Matrix aux_X {config.nodes, std::vector<double>(config.nodes, 0.0)};
    Matrix aux_Y {config.nodes, std::vector<double>(config.nodes, 0.0)};

    Matrix gradPhi_X {config.nodes, std::vector<double>(config.nodes, 0.0)};
    Matrix gradPhi_Y {config.nodes, std::vector<double>(config.nodes, 0.0)};

    Matrix divAuxField {config.nodes, std::vector<double>(config.nodes, 0.0)};

    Matrix dux_dx {config.nodes, std::vector<double>(config.nodes, 0.0)};
    Matrix dux_dy {config.nodes, std::vector<double>(config.nodes, 0.0)};

    Matrix duy_dx {config.nodes, std::vector<double>(config.nodes, 0.0)};
    Matrix duy_dy {config.nodes, std::vector<double>(config.nodes, 0.0)};

    Matrix dd_ux {config.nodes, std::vector<double>(config.nodes, 0.0)};
    Matrix dd_uy {config.nodes, std::vector<double>(config.nodes, 0.0)};

    Matrix lap_phi {config.nodes, std::vector<double>(config.nodes, 0.0)};

    while (counter < config.endTime)
    {

        for (int i {0} ; i < (config.nodes) ; ++i)
        {
            aux_X[0][i] = botx[i];
            aux_Y[0][i] = boty[i];
            aux_X[config.nodes - 1][i] = upx[i];
            aux_Y[config.nodes - 1][i] = upy[i];
            aux_X[i][0] = leftx[i];
            aux_Y[i][0] = lefty[i];
            aux_X[i][config.nodes - 1] = rightx[i];
            aux_Y[i][config.nodes - 1] = righty[i];
        }

        for (int i {1} ; i < (config.nodes - 2) ; ++i)
        {
            for (int j {1} ; j < (config.nodes - 2) ; ++j)
            {
                dux_dx[i][j] = (u_x[i][j + 1] - u_x[i][j - 1]) / (2.0 * config.dx);
                dux_dy[i][j] = (u_x[i + 1][j] - u_x[i - 1][j]) / (2.0 * config.dx); // dx == dy so dx and dy interchangeable

                duy_dx[i][j] = (u_y[i][j + 1] - u_y[i][j - 1]) / (2.0 * config.dx);
                duy_dy[i][j] = (u_y[i + 1][j] - u_y[i - 1][j]) / (2.0 * config.dx); 

                dd_ux[i][j] = ((u_x[i][j + 1] - 2.0 * u_x[i][j] + u_x[i][j - 1]) / (config.dx * config.dx)) + ((u_x[i + 1][j] - 2.0 * u_x[i][j] + u_x[i - 1][j]) / (config.dx * config.dx));
                dd_uy[i][j] = ((u_y[i][j + 1] - 2.0 * u_y[i][j] + u_y[i][j - 1]) / (config.dx * config.dx)) + ((u_y[i + 1][j] - 2.0 * u_y[i][j] + u_y[i - 1][j]) / (config.dx * config.dx));

                aux_X[i][j] = u_x[i][j] + dt * (-1.0 * (u_x[i][j] * dux_dx[i][j] + u_y[i][j] * dux_dy[i][j]) + config.visc * dd_ux[i][j]);
                aux_Y[i][j] = u_y[i][j] + dt * (-1.0 * (u_x[i][j] * duy_dx[i][j] + u_y[i][j] * duy_dy[i][j]) + config.visc * dd_uy[i][j]);
            
                divAuxField[i][j] = ((aux_X[i][j + 1] - aux_X[i][j]) / config.dx + (aux_Y[i + 1][j] - aux_Y[i][j]) / config.dx);
            }
        }
        

        // Poisson solve
        double phi_counter {0.0};
        double residual_max {std::numeric_limits<double>::infinity()};

        while (1e-5 < residual_max && phi_counter < static_cast<double>(config.phi_iterations))
        {
            Matrix old_phi {phi};

            for (int i {1} ; i < (config.nodes - 2) ; ++i)
            {
                for (int j {1} ; j < (config.nodes - 1) ; ++j)
                {
                    phi[i][j] = (0.25) * (old_phi[i + 1][j] + old_phi[i - 1][j] + old_phi[i][j + 1] + old_phi[i][j - 1] - ((config.dx * config.dx) * divAuxField[i][j]));
                    
                    lap_phi[i][j] = ((phi[i][j + 1] - 2 * phi[i][j] + phi[i][j - 1]) / (config.dx * config.dx) + (phi[i + 1][j] - 2 * phi[i][j] + phi[i - 1][j]) / (config.dx * config.dx));
                }
            }

            residual_max = maxMatrix(absMatrix(subtractMatrices(lap_phi, divAuxField)));
            phi_counter += 1;
        }

        // Calculate gradient of phi using backward differences
        for (int i {1} ; i < (config.nodes - 1) ; ++i)
        {
            for (int j {1} ; j < (config.nodes - 1) ; ++j)
            {
                gradPhi_X[i][j] = (phi[i][j] - phi[i][j - 1]) / config.dx;
                gradPhi_Y[i][j] = (phi[i][j] - phi[i - 1][j]) / config.dx;
            }
        }

        // Projection step
        u_x = subtractMatrices(aux_X, gradPhi_X);
        u_y = subtractMatrices(aux_Y, gradPhi_Y);

        // Reapply boundary conditions
        for (int i {0} ; i < config.nodes ; ++i)
        {
            u_x[0][i] = botx[i];
            u_y[0][i] = boty[i];
            u_x[config.nodes - 1][i] = upx[i];
            u_y[config.nodes - 1][i] = upy[i];
            u_x[i][0] = leftx[i];
            u_y[i][0] = lefty[i];
            u_x[i][config.nodes - 1] = rightx[i];
            u_y[i][config.nodes - 1] = righty[i];
        }

        counter += dt;

        double speed_x {maxMatrix(absMatrix(u_x))};
        double speed_y {maxMatrix(absMatrix(u_y))};

        if (speed_x == 0 && speed_y == 0)
        {
            advective = std::numeric_limits<double>::infinity();
        }
        else
        {
            advective = config.cfl / ((speed_x / config.dx) + (speed_y / config.dx));
        }

        diffusive = (config.dx * config.dx) / (4 * config.visc);
        dt = std::min(advective, diffusive);

        std::cout << "dt: " << dt << "\ndiffusive: " << diffusive << "\n";
    }
    return 0;
}
