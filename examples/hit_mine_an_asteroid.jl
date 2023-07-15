using GTOC12
using Plots
using LinearAlgebra
using Interpolations

# initialize()
# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1))
asteroid = Asteroid(ID_min)
furnish_all_kernels()
DV₀ = [-1000; 0; 3500]
ET_start = GTOC12.ET₀ + 340*24*3600
x_Earth = get_body_state(GTOC12.Earth; ET=ET_start)
x₀ = x_Earth + [0; 0; 0; DV₀[:]]

# Transfer time
Δt = 1/3 * 365 * 24 * 3600.0
ΔV∞_max = 6000.0;

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = ET_start + Δt
x_target = get_body_state(asteroid; ET=ET_target)
r_target = view(x_target, 1:3)

# NOTE: This is super brittle to initial conditions
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
# v⃗₀, v⃗ = lambert(; r⃗₀=x₀[1:3], r⃗=r_target, Δt, tₘ=1)
# x₀⁺ = vcat(x₀[1:3], v⃗₀)
# xₜ = vcat(r_target,v⃗) 
x₀⁺, xₜ = fixed_time_single_shoot(x₀, Δt, r_target; print_iter=true)
ΔV_departure_1 = (x₀⁺-x_Earth)[4:6]
println("ΔV_departure_1 $(norm(ΔV_departure_1))")
ΔV_arrival_1 = (xₜ-x_target)[4:6]

# Initialize Mining Ship
mining_ship = GTOC12.Mining_Ship()

line_array = Vector{String}()
line_array, mining_ship = GTOC12.record_line(line_array, "launch", [x_Earth, x₀⁺], ET_start, mining_ship)

# hit an asteroid 
t0 = ET_start
every_day = 24 * 60 * 60
x_spacecraft, T_spacecraft, time_ET = calculate_rendezvous_from_Earth(x₀⁺, x_target, Δt, t0; m0=mining_ship.mass_total,
                                                                    μ=GTOC12.μ_☉, dt=24*3600, output_times=every_day,
                                                                    abstol = 1e-4, reltol=1e-6)
r_burn_arc_1 = getindex.(x_spacecraft', 1:3)'
println("asteroid_rendezvous_resid $(x_spacecraft[end] - x_target)")


line_array, mining_ship = GTOC12.record_line(line_array, "burn", x_spacecraft, time_ET, mining_ship, control=T_spacecraft)
# Deploy a miner
# TODO need to update mass in the mining_ship so that it can be used in this file 
# just output mining structure as well 
# could change this to event_line? 
# TODO mass is wrong in the second line when deploying 
# 
# hit an asteroid, deploy a miner
line_array, mining_ship = GTOC12.record_line(line_array, "rendezvous", x_spacecraft[end], time_ET[end], mining_ship, rendez_flag="deploy", event_ID=ID_min)

# hit an asteroid, recover a miner
tf = time_ET[end] + every_day*10
line_array, mining_ship = GTOC12.record_line(line_array, "rendezvous", x_spacecraft[end], tf, mining_ship, rendez_flag="recover", event_ID=ID_min)

# # rendezvous with earth 




# write lines to txt file
println("Printing results to file.")
file = open("Result.txt", "w")
for line in line_array
    write(file, line)
end
close(file)

# # Remove empty lines from file
# run(`sed Result.txt`)

# # Move file into problem directory
# mv("Result.txt", "../problem/GTOC12_Verification_Program/GTOC12_Verification/Linux/")

# Run verification program
# run(``)

# Plotting
# Plot Earth/asteroid/spacecraft
ETs = [ET_start, ET_target]
plot_body(asteroid; ETs=ETs, color="Red")
plot_body!(GTOC12.Earth; ETs=ETs, color="Blue")
plot_coast!(x₀⁺, Δt; label="Coast 1", color="Green")
plot!(r_burn_arc_1[:,1], r_burn_arc_1[:,2], r_burn_arc_1[:,3], color="cyan", label="steered burn 1")
