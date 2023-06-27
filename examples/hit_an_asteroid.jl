using GTOC12
using Plots

# Import asteroids
asteroid_df = get_asteroid_df()

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(asteroid_df.sma/au2m .- 1))
asteroid = asteroid_df[ID_min, :]

# Try out shooting method to hit asteroid
DV₀ = [-1000; 0; 0]
x₀ = get_planet_state("EARTH", GTOC12.ET₀) + [0;0;0;DV₀[:]]

# Transfer time
t_transfer = 1/2 * 365 * 24 * 3600

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = GTOC12.ET₀ + t_transfer
x_target = get_asteroid_state(asteroid, ET_target)
r_target = view(x_target, 1:3)

# NOTE: This is super brittle to initial conditions
# TODO: add graceful error handling
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
x₀⁺, xₜ = fixed_time_single_shoot( x₀, t_transfer, r_target; print_iter=false)
ΔV_departure_1 = ( x₀⁺ - x₀ )[4:6]
ΔV_arrival_1 = (xₜ - x_target )[4:6]

# Add small coast period
t_coast = 24*3600

# Try to hit Earth again
x₀ = xₜ
t_transfer2 = t_transfer

ET_arrival = GTOC12.ET₀ + t_transfer + t_coast + t_transfer2
x_target = get_planet_state("EARTH", ET_arrival)
r_target = x_target[1:3]
x₀⁺_2, xₜ = fixed_time_single_shoot( x₀, t_transfer, r_target; print_iter=false)
ΔV_departure_2 = ( x₀⁺_2 - x₀ )[4:6]
ΔV_arrival_2 = (xₜ - x_target )[4:6]



# Plot Earth/asteroid/spacecraft
ETs = [GTOC12.ET₀, ET_arrival]
plot_asteroid_from_df_row(asteroid; ET_in=ETs)
plot_planet!(planet="EARTH"; ET_in=ETs, label="Earth", color=colormap("Blues"))
plot_coast!(x₀⁺, t_transfer; label="Coast 1", color=colormap("Greens"))
plot_coast!(x₀⁺_2, t_transfer2; label="Coast 2", color=colormap("Greens"))

# Holy shit, this is actually a feasible solution we can use for scoring
# Let's try to write it to a scorable file
