using GTOC12
using Plots

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1))
asteroid = Asteroid(ID_min)

# Try out shooting method to hit asteroid
furnish_all_kernels()
DV₀ = [-1000; 0; 0]
x₀ = get_body_state(GTOC12.Earth; ET=GTOC12.ET₀) + [0; 0; 0; DV₀[:]]

# Transfer time
Δt = 1 / 2 * 365 * 24 * 3600
ΔV∞_max = 6000.0;

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = GTOC12.ET₀ + Δt
x_target = get_body_state(asteroid; ET=ET_target)
r_target = view(x_target, 1:3)

# NOTE: This is super brittle to initial conditions
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
x₀⁺, xₜ = fixed_time_single_shoot(x₀, Δt, r_target; print_iter=false)
ΔV_departure_1 = (x₀⁺-x₀)[4:6]
ΔV_arrival_1 = (xₜ-x_target)[4:6]

# Plot Earth/asteroid/spacecraft
ETs = [GTOC12.ET₀, ET_target]
# plot_asteroid_from_df_row(asteroid; ET_in=ETs)
plot_body(asteroid; ETs=ETs, color=colormap("Reds"))
plot_body!(GTOC12.Earth; ETs=ETs, color=colormap("Blues"))
plot_coast!(x₀⁺, Δt; label="Coast 1", color=colormap("Greens"))

# TODO: Optimize the above for DV/via Time-of-flight

# TODO: For launch case, just need initial time, not states
# states_out, controls_out, time_out = optimize_impulsive_launch(x₀, x_target, Δt; ΔV₀=ΔV_departure_1, asteroid_ID=ID_min)


# Try running continuous burn over this arc
states_out, controls_out, time_out = optimize_continuous_arc(x₀⁺, x_target, Δt; asteroid_ID=ID_min, m₀=500.0)

