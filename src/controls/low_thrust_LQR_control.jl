using GTOC12
using SPICE
using ForwardDiff
using LinearAlgebra
using Plots

# Import asteroids
asteroid_df = get_asteroid_df()
using Plots

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1))
asteroid = GTOC12.asteroid_df[ID_min, :]

# Try out shooting method to hit asteroid
furnish_all_kernels()
DV₀ = [-1000; 0; 0]
x₀ = get_planet_state("EARTH", GTOC12.ET₀) + [0; 0; 0; DV₀[:]]

# Transfer time
Δt = 1 / 2 * 365 * 24 * 3600
Δt = 24 * 3600
ΔV∞_max = 6000.0;

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = GTOC12.ET₀ + Δt
x_target = get_asteroid_state(asteroid, ET_target)
println("x_target: ")
println(x_target)
r_target = view(x_target, 1:3)

# NOTE: This is super brittle to initial conditions
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
x₀⁺, xₜ = fixed_time_single_shoot(x₀, Δt, r_target; print_iter=false)
ΔV_departure_1 = (x₀⁺-x₀)[4:6]
ΔV_arrival_1 = (xₜ-x_target)[4:6]

# Add small coast period
t_coast = 24 * 3600

# # Try to hit Earth again
# x₀ = xₜ
# t_transfer2 = t_transfer

# ET_arrival = GTOC12.ET₀ + t_transfer + t_coast + t_transfer2
# x_target = get_planet_state("EARTH", ET_arrival)
# r_target = x_target[1:3]
# x₀⁺_2, xₜ = fixed_time_single_shoot( x₀, t_transfer, r_target; print_iter=false)
# ΔV_departure_2 = ( x₀⁺_2 - x₀ )[4:6]
# ΔV_arrival_2 = (xₜ - x_target )[4:6]


# Plot Earth/asteroid/spacecraft
ETs = [GTOC12.ET₀, ET_target]
plot_asteroid_from_df_row(asteroid; ET_in=ETs)
plot_planet!(planet="EARTH"; ET_in=ETs, label="Earth", color=colormap("Blues"))
plot_coast!(x₀⁺, Δt; label="Coast 1", color=colormap("Greens"))
# plot_coast!(x₀⁺_2, t_transfer2; label="Coast 2", color=colormap("Greens"))

# TODO: Optimize the above for DV
# TODO: Create an abstract type for both planets and asteroids
# Since they're both on Keplerian rails, but planets can conduct flybys and asteroids can provide resources

# TODO: For launch case, just need initial time, not states
# update dis????
states_out, controls_out, time_out = optimize_continuous_launch(x₀, x_target, Δt; ΔV₀=ΔV_departure_1, asteroid_ID=ID_min)

