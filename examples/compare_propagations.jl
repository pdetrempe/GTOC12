using GTOC12
using Plots
using TrajectoryOptimization
using RobotDynamics
using RobotDynamics: @autodiff
using StaticArrays
using ForwardDiff, FiniteDiff
using LinearAlgebra

const TO = TrajectoryOptimization
const RD = RobotDynamics

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma/au2m .- 1))
asteroid = GTOC12.asteroid_df[ID_min, :]

# Try out shooting method to hit asteroid
furnish_all_kernels()
DV₀ = [0; 0; 0]
x₀ = get_planet_state("EARTH", GTOC12.ET₀) + [0;0;0;DV₀[:]]

# Transfer time
t_transfer = 1/2 * 365 * 24 * 3600

# Transfer via Keplerian propagation
x_transfer_Kepler = propagate_keplerian(x₀, t_transfer)
print(x₀)

# Transfer via universal propagation
x_transfer_universal= propagate_universal(x₀, t_transfer)
print(x₀)


# Plot Earth/asteroid/spacecraft
ETs = [GTOC12.ET₀, ET_target]
plot_asteroid_from_df_row(asteroid; ET_in=ETs)
plot_planet!(planet="EARTH"; ET_in=ETs, label="Earth", color=colormap("Blues"))
plot_coast!(x_transfer_Kepler, t_transfer; label="Kepler", color="cyan")
plot_coast!(x_transfer_universal, t_transfer; label="Universal", color="magenta")