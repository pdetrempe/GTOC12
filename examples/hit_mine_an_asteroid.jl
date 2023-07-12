using GTOC12
using Plots
using LinearAlgebra

# initialize()
# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1))
asteroid = Asteroid(ID_min)
furnish_all_kernels()
DV₀ = [-1000; 0; 3500]
x₀ = get_body_state(GTOC12.Earth; ET=GTOC12.ET₀) + [0; 0; 0; DV₀[:]]

# Transfer time
Δt = 3/4 * 365 * 24 * 3600
ΔV∞_max = 6000.0;

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = GTOC12.ET₀ + Δt
x_target = get_body_state(asteroid; ET=ET_target)
r_target = view(x_target, 1:3)

# NOTE: This is super brittle to initial conditions
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
# v⃗₀, v⃗ = lambert(; r⃗₀=x₀[1:3], r⃗=r_target, Δt, tₘ=1)
# x₀⁺ = vcat(x₀[1:3], v⃗₀)
# xₜ = vcat(r_target,v⃗) 
x₀⁺, xₜ = fixed_time_single_shoot(x₀, Δt, r_target; print_iter=true)
ΔV_departure_1 = (x₀⁺-x₀)[4:6]
ΔV_arrival_1 = (xₜ-x_target)[4:6]

# Initialize Mining Ship
mining_ship = GTOC12.Mining_Ship()

line_array = GTOC12.record_line("launch", x₀, x₀⁺, GTOC12.ET₀, mining_ship, 0.0)


# hit an asteroid 

# deploy a miner 

# recover a miner

# rendezvous with earth 

# output dataframe

