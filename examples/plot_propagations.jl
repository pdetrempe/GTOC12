using GTOC12
using LinearAlgebra
using Plots

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma/au2m .- 1))
asteroid = GTOC12.asteroid_df[ID_min, :]

ET_start = GTOC12.ET₀ + 1/4*365*24*3600

# Try out shooting method to hit asteroid
furnish_all_kernels()
DV₀ = [0; 800; 0]
x₀ = get_planet_state("EARTH", ET_start) + [0;0;0;DV₀[:]]

# Transfer time
t_transfer = 1/3 * 365 * 24 * 3600
ET_target = ET_start + t_transfer
x_target = get_asteroid_state(asteroid, ET_target)
r_target = view(x_target, 1:3)

# Try
# Transfer via universal propagation
x_transfer_universal= propagate_universal(x₀, t_transfer)

# NOTE: This is super brittle to initial conditions
x₀⁺, xₜ = fixed_time_single_shoot( x₀, t_transfer, r_target; print_iter=true)
ΔV_departure_shooting = ( x₀⁺ - x₀ )[4:6]
ΔV_arrival_shooting = (xₜ - x_target )[4:6]

# TODO: do a time sweep of Lambert trajectories to minimize DV
# Lambert solution
v⃗₀, v⃗ = lambert(; r⃗₀=x₀[1:3], r⃗=r_target, Δt=t_transfer, tₘ=1) 
x_lambert = vcat(x₀[1:3], v⃗₀)
ΔV_departure_lambert = v⃗₀ - x₀[4:6]
ΔV_arrival_lambert = v⃗ - x_target[4:6]

# Plot Earth/asteroid/spacecraft
ETs = [GTOC12.ET₀, ET_target]
plot_asteroid_from_df_row(asteroid; ET_in=ETs)
# plot_planet!(planet="EARTH"; ET_in=ETs, label="Earth", color=colormap("Blues"))
plot_coast!(x₀⁺, t_transfer; label="Shooting Method", color="cyan")
plot_coast!(x_lambert, t_transfer; label="Lambert Method", color="black")
plot_coast!(x₀, t_transfer; label="Propagation", color="magenta")