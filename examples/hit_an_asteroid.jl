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
r_target = get_asteroid_state(asteroid, ET_target)[1:3]

# NOTE: This is super brittle to initial conditions
# TODO: add graceful error handling
x₀⁺ = fixed_time_single_shoot( x₀, t_transfer, r_target; print_iter=false)


# Plot Earth/asteroid/spacecraft
ETs = [GTOC12.ET₀, ET_target]
plot_asteroid_from_df_row(asteroid; ET_in=ETs)
plot_planet!(planet="EARTH"; ET_in=ETs, label="Earth", color=colormap("Blues"))
plot_coast!(x₀⁺, t_transfer; label="coast", color=colormap("Greens"))

