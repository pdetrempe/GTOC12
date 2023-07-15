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

line_array = []
line_array, mining_ship = GTOC12.record_line(line_array, "launch", [x₀, x₀⁺], GTOC12.ET₀, mining_ship)

# hit an asteroid 
t0 = GTOC12.ET₀
every_day = 24 * 60 * 60
x_spacecraft, T_spacecraft, time_ET = calculate_rendezvous(x₀⁺, x_target, Δt, t0; m0=mining_ship.mass_total,
                                                                    μ=GTOC12.μ_☉, dt=24*3600, output_times=every_day)


# x_burn_arc_1, T_vector, time_vector_ET = calculate_rendezvous(x₀⁺, x_target, Δt, GTOC12.ET₀; m0=m, μ=GTOC12.μ_☉, dt=24*3600,
# abstol=1e-10, 
# reltol=1e-14)
#for (c,d) in zip(time_ET, eachrow(T_spacecraft))
    #println(d)
    #println(length(d))
#end
#x_spacecraft_updated = interpolate(time_ET, x_spacecraft, desired_time_ET, Gridded(Linear()))

line_array, mining_ship = GTOC12.record_line(line_array, "burn", x_spacecraft, time_ET, mining_ship, control=T_spacecraft)
# Deploy a miner
# TODO need to update mass in the mining_ship so that it can be used in this file 
# just output mining structure as well 
# could change this to event_line? 
# TODO mass is wrong in the second line when deploying 

# hit an asteroid, deploy a miner
line_array, mining_ship = GTOC12.record_line(line_array, "rendezvous", x_spacecraft[end], time_ET[end], mining_ship, rendez_flag="deploy", event_ID=ID_min)


# recover a miner

# rendezvous with earth 




# write lines to txt file
file = open("Result.txt", "w")
for line in line_array
    write(file, line)
end
close(file)