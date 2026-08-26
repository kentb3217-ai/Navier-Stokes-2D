import matplotlib.pyplot as plt, numpy as np

# Defining the equation
length = 1 # Keep length at 1, anything greater than 3 and it gets shaky
nodes = 20
time = 2
C = 0.1 # CFL number, typically between 0 - 1
velocity_bound_cond = 1 # ideally between 0 - 1
visc = .01 # viscocity, typically between 0.01 - .1, in more complex systems it isn't a constant (depends on temp), .01 is what i normally put it at
phi_iterations = 500 # number of iterations phi can be calculated and refined

dx = length / (nodes - 1)
dy = length / (nodes - 1)

# velocity of fluid in x and y directions respectively
u_x = np.zeros((nodes, nodes)) # u
u_y = np.zeros((nodes, nodes)) # v in many cases

# boundary conditions
botx = 0
boty = 0
upx = np.full(nodes, velocity_bound_cond)
upy = 0
leftx = 0
lefty = 0
rightx = 0
righty = 0

u_x[0, :] = botx
u_y[0, :] = boty
u_x[-1, :] = upx
u_y[-1, :] = upy
u_x[:, 0] = leftx
u_y[:, 0] = lefty
u_x[:, -1] = rightx
u_y[:, -1] = righty

# Time step
# advective restriction (bad if max velo is zero), hence why we add the diffusive restriction
speed_x = np.max(np.abs(u_x))
speed_y = np.max(np.abs(u_y))

if (speed_x == 0 and speed_y == 0):
    dt_adv = np.inf
else:
    dt_adv = C / ((speed_x / dx) + (speed_y / dy))

# diffusive restriction
dt_diff = dx**2 / (4 * visc)
dt = min(dt_adv, dt_diff)

# Initialization for phi and aux_field_x and aux_field_y
phi = np.zeros((nodes, nodes))
aux_field_x = np.zeros((nodes, nodes))
aux_field_y = np.zeros((nodes, nodes))

aux_field_x[0, :] = botx
aux_field_y[0, :] = boty
aux_field_x[-1, :] = upx
aux_field_y[-1, :] = upy
aux_field_x[:, 0] = leftx
aux_field_y[:, 0] = lefty
aux_field_x[:, -1] = rightx
aux_field_y[:, -1] = righty

grad_phi_x = np.zeros((nodes, nodes))
grad_phi_y = np.zeros((nodes, nodes))
div_aux_field = np.zeros((nodes, nodes))

# Counter
counter = 0

# for plotting
fig, ax = plt.subplots()
x = np.arange(nodes)
y = np.arange(nodes)
X, Y = np.meshgrid(x, y)

while counter < time:
    
    # Calculate auxiliary field of u (velocity)
    dux_dx = (u_x[1:-1, 2:] - u_x[1:-1, :-2]) / (2*dx)
    dux_dy = (u_x[2:, 1:-1] - u_x[:-2, 1:-1]) / (2*dy)

    duy_dx = (u_y[1:-1, 2:] - u_y[1:-1, :-2]) / (2*dx)
    duy_dy = (u_y[2:, 1:-1] - u_y[:-2, 1:-1]) / (2*dy)

    # Calculate the Laplacian of u
    dd_ux = ((u_x[1:-1, 2:] - 2 * u_x[1:-1, 1:-1] + u_x[1:-1, :-2]) / dx**2) + ((u_x[2:, 1:-1] - 2 * u_x[1:-1, 1:-1] + u_x[:-2, 1:-1]) / dy**2)
    dd_uy = ((u_y[1:-1, 2:] - 2 * u_y[1:-1, 1:-1] + u_y[1:-1, :-2]) / dx**2) + ((u_y[2:, 1:-1] - 2 * u_y[1:-1, 1:-1] + u_y[:-2, 1:-1]) / dy**2)

    # Calculate auxiliary field
    aux_field_x[1:-1, 1:-1] = u_x[1:-1, 1:-1] + dt * (-1 * (u_x[1:-1, 1:-1] * dux_dx + u_y[1:-1, 1:-1] * dux_dy) + visc * dd_ux)
    aux_field_y[1:-1, 1:-1] = u_y[1:-1, 1:-1] + dt * (-1 * (u_x[1:-1, 1:-1] * duy_dx + u_y[1:-1, 1:-1] * duy_dy) + visc * dd_uy)

    aux_field_x[0, :] = botx
    aux_field_y[0, :] = boty
    aux_field_x[-1, :] = upx
    aux_field_y[-1, :] = upy
    aux_field_x[:, 0] = leftx
    aux_field_y[:, 0] = lefty
    aux_field_x[:, -1] = rightx
    aux_field_y[:, -1] = righty

    # Calculate the divergence of the auxiliary field
    div_aux_field[1:-1, 1:-1] = ((aux_field_x[1:-1, 2:] - aux_field_x[1:-1, 1:-1]) / dx + (aux_field_y[2:, 1:-1] - aux_field_y[1:-1, 1:-1]) / dy)

    # Poisson solve
    phi_counter = 0
    residual_max = np.inf
    while 1e-5 < residual_max and phi_counter < phi_iterations:
        # Store old phi
        old_phi = np.copy(phi)

        # Calculate phi utilizing Jacobi iteration
        phi[1:-1, 1:-1] = (0.25) * (old_phi[2:, 1:-1] + old_phi[:-2, 1:-1] + old_phi[1:-1, 2:] + old_phi[1:-1, :-2] - (dx**2 * div_aux_field[1:-1, 1:-1]))

        # boundary conditions for phi (dependent on inside values) (no change in phi)
        phi[0, :] = phi[1, :]
        phi[-1, :] = phi[-2, :]
        phi[:, 0] = phi[:, 1]
        phi[:, -1] = phi[:, -2]
        phi[0, 0] = 0 # fix value, since phi is only defined up to a constant

        # Calculate residual
        lap_phi = ((phi[1:-1, 2:] - 2 * phi[1:-1, 1:-1] + phi[1:-1, :-2]) / dx**2 + (phi[2:, 1:-1] - 2 * phi[1:-1, 1:-1] + phi[:-2, 1:-1]) / dy**2)
        residual_max = np.max(np.abs(lap_phi - div_aux_field[1:-1, 1:-1]))

        phi_counter += 1

    # Calculate gradient of phi using backward differences
    grad_phi_x[1:-1, 1:-1] = (phi[1:-1, 1:-1] - phi[1:-1, :-2]) / dx
    grad_phi_y[1:-1, 1:-1] = (phi[1:-1, 1:-1] - phi[:-2, 1:-1]) / dy

    # Projection step
    u_x = aux_field_x - grad_phi_x
    u_y = aux_field_y - grad_phi_y

    # reapply boundary conditions
    u_x[0, :] = botx
    u_y[0, :] = boty
    u_x[-1, :] = upx
    u_y[-1, :] = upy
    u_x[:, 0] = leftx
    u_y[:, 0] = lefty
    u_x[:, -1] = rightx
    u_y[:, -1] = righty

    ax.clear()
    ax.set_aspect('equal')
    ax.quiver(X, Y, u_x, u_y)
    ax.set_title("Distribution at t: {:.3f} [s].".format(counter))
    plt.pause(0.01)

    counter += dt

    # calculate max velocity in x and y directions
    speed_x = np.max(np.abs(u_x))
    speed_y = np.max(np.abs(u_y))

    # Calculate timestep
    if (speed_x == 0 and speed_y == 0):
        dt_adv = np.inf
    else:
        dt_adv = C / ((speed_x / dx) + (speed_y / dy))
    dt_diff = dx**2 / (4 * visc)
    dt = min(dt_adv, dt_diff)

    # Final divergence diagnostic
    div_u = ((u_x[1:-1, 2:] - u_x[1:-1, 1:-1]) / dx + (u_y[2:, 1:-1] - u_y[1:-1, 1:-1]) / dy)

    # For debugging purposes
    print(f"""dt: {dt},
        u_x: {np.max(np.abs(u_x))},
        u_y: {np.max(np.abs(u_y))},
        div_aux: {np.max(np.abs(div_aux_field))}
        aux_fieldx: {np.max(np.abs(aux_field_x))}
        aux_fieldy: {np.max(np.abs(aux_field_y))}
        projected_div: {np.max(np.abs(div_u))}""")

plt.show()
