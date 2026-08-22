import matplotlib.pyplot as plt, numpy as np

# Defining the equation
length = 1 # Keep length at 1, anything greater than 3 and it gets shaky
nodes = 20
time = 100
C = 0.1 # CFL number, typically between 0 - 1
velocity_bound_cond = 1 # ideally between 0 - 1
visc = .01 # viscocity, typically between 0.01 - .1, in more complex systems it isn't a constant (depends on temp), .01 is what i normally put it at
phi_iterations = 500 # number of iterations phi can be calculated and refined
velocity_threshold = velocity_bound_cond * 10**-6 # To ensure velocity doesn't go below a certain value, otherwise sim will blow up

dx = length / (nodes - 1)
dy = length / (nodes - 1)

# velocity of fluid in x and y directions respectively
u_x = np.zeros((nodes, nodes)) # u
u_y = np.zeros((nodes, nodes)) # v in many cases

# boundary conditions
botx = np.linspace(0, -velocity_bound_cond, nodes)
boty = np.linspace(0, velocity_bound_cond, nodes)
upx = 0
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
dt_adv = C * min(dx/(np.max(np.abs(u_x))), dy/(np.max(np.abs(u_y)))) 
# take np.max(np.abs(u)) b/c we want to know the fastest something is moving and move time_step according to that
# otherwise, if we don't, then we end up with messy garbage because time step didn't take into account the fastest moving point in the navier stokes

# diffusive restriction
dt_diff = dx**2 / (4 * visc)
dt = min(dt_adv, dt_diff)

# Initialization for phi and aux_field_x and aux_field_y
phi = np.zeros((nodes, nodes))
aux_field_x = np.zeros((nodes, nodes))
aux_field_y = np.zeros((nodes, nodes))
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
    phi_counter = 0
    change_phi = np.inf
    
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

    # Calculate the divergence of the auxiliary field
    div_aux_field[1:-1, 1:-1] = ((aux_field_x[1:-1, 2:] - aux_field_x[1:-1, :-2]) / (2*dx)) + ((aux_field_y[2:, 1:-1] - aux_field_y[:-2, 1:-1]) / (2*dy))
    
    # Poisson solve
    while 1e-5 < change_phi and phi_counter < phi_iterations:
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

        # Calculate change in phi, max change 10 ** -5
        change_phi = np.max(np.abs(phi - old_phi))

        phi_counter += 1

    # Calculate gradient of phi, utilized central first differences formula
    grad_phi_x[1:-1, 1:-1] = (phi[1:-1, 2:] - phi[1:-1, :-2]) / (2*dx)
    grad_phi_y[1:-1, 1:-1] = (phi[2:, 1:-1] - phi[:-2, 1:-1]) / (2*dy)

    # Find grad of phi at boundaries using forward and backward differences 
    # (top and right boundary --> backward, bottom and left boundary --> forward)
    # top
    grad_phi_y[0, :] = (phi[1, :] - phi[0, :]) / dy 
    # bottom
    grad_phi_y[-1, :] = (phi[-1, :] - phi[-2, :]) / dy
    # left
    grad_phi_x[:, 0] = (phi[:, 1] - phi[:, 0]) / dx
    # right
    grad_phi_x[:, -1] = (phi[:, -1] - phi[:, -2]) / dx

    # compute velocity vector for u_x and u_y for n+1
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
    vx_max = np.max(np.abs(u_x))
    vy_max = np.max(np.abs(u_y))

    # Calculate timestep
    dt_adv = C * min( dx / vx_max, dy / vy_max )
    dt_diff = dx**2 / (4 * visc)
    dt = min(dt_adv, dt_diff)

    # For debugging purposes
    print(f"""dt: {dt},
        u_x: {np.max(np.abs(u_x))},
        u_y: {np.max(np.abs(u_y))},
        div_aux: {np.max(np.abs(div_aux_field))}
        aux_fieldx: {np.max(np.abs(aux_field_x))}
        aux_fieldy: {np.max(np.abs(aux_field_y))}""")

plt.show()
