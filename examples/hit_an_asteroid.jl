using GTOC12
using SPICE
using ForwardDiff
using LinearAlgebra

# Import asteroids
asteroid_df = get_asteroid_df()

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(asteroid_df.sma .- 1))

# Plot the asteroid trajectory
asteroid = asteroid_df[ID_min, :]
plot_asteroid_from_df_row(asteroid)

# Plot Earth
plot_planet!(planet="EARTH"; label="Earth", color=colormap("Blues"))

# Try out shooting method to hit asteroid
DV₀ = [-GTOC12.v∞_max; 0; 0]
x₀ = get_planet_state("EARTH", GTOC12.ET₀) + [0;0;0;DV₀[:]]

t_shoot = 1 / 2 * 365 * 24 * 3600

# Fixed-time shooting
# 1. Get asteroid position (target)
ET_target = GTOC12.ET₀ + t_shoot
x_target = get_asteroid_state(asteroid, ET_target)

# 2. Propagate spacecraft from initial state
xₜ = propagate_keplerian(x₀, t_shoot)

# 3. Get STM: use ForwardDiff to get dr_f/dv_0
# Actual "cost" is miss distance
# Need f(v₀) = dr
function calculate_asteroid_miss_distance(x_aug)
    # Initial guess state values
    x₀ = x_aug[1:6]
    t = x_aug[7]
    rₜ_des = x_aug[8:10]

    # TODO: Augment target state to calculate sensitivity w.r.t. initial ET

    # ET_target = GTOC12.ET₀ + t
    # asteroid = get_asteroid_df()[ID, :]
    # x_target = get_asteroid_state(asteroid, ET_target)
    xₜ = propagate_keplerian(x₀, t)

    cost = xₜ[1:3] - rₜ_des
end

# TODO: Think about clever way to package problem into a type
# such that can have one giant vector for calculating adjoint_sensitivities
# 
# 5. Iterate until convergence
MAX_ITER = 50
num_iter = 0
drₜ = x₀[1:3]
while num_iter < MAX_ITER && norm(drₜ) > GTOC12.POS_ABS_TOL
    # Form vector used by ForwardDiff problem
    x_aug = vcat(x₀, t_shoot, x_target[1:3])

    # Calculate miss distance from current propagation
    global drₜ = calculate_asteroid_miss_distance(x_aug)
    println(drₜ)

    # Calculate sensitivity of solution to all parameters
    dr_dx_aug = ForwardDiff.jacobian(calculate_asteroid_miss_distance, x_aug)

    # Get partials w.r.t. design variables
    # (In this case, initial velocity and time of flight)
    ∂rₜ_∂v₀ = dr_dx_aug[:, 4:6] # Variation of final position w.r.t. initial velocity
    ∂rₜ_∂TOF = dr_dx_aug[:, 7]  # Variation of final position w.r.t. time-of-flight

    # ^^^TODO: handle indices more intelligently for packing/unpacking


    # Update initial state via Newton's method
    # See Pavlak Thesis Eq. 3.12
    dv₀ = -∂rₜ_∂v₀\drₜ
    dTOF = -∂rₜ_∂TOF\drₜ

    global x₀[4:6] += dv₀

    # Variable time stuff isn't working interestingly
    # global t_shoot += dTOF

    global num_iter += 1

end

plot_coast!(x₀, t_shoot; label="coast", color=colormap("Greens"))
