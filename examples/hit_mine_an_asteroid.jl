using GTOC12
using Plots
using LinearAlgebra
using Interpolations
using OrdinaryDiffEq

# Initialize Mining Ship
mining_ship = GTOC12.Mining_Ship()

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1))
asteroid = Asteroid(ID_min)
furnish_all_kernels()
DV₀ = [159.31762519335462; 1505.6951917938995; 3569.929588119754]
ET_start = GTOC12.ET₀ + 338*24*3600
x_Earth = get_body_state(GTOC12.Earth; ET=ET_start)
x₀ = x_Earth + [0; 0; 0; DV₀[:]]

# Transfer time
# Δt = 0.3207 * 365 * 24 * 3600.0
Δt = 0.35 * 365*24*3600
ΔV∞_max = 6000.0;

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = ET_start + Δt
x_target = get_body_state(asteroid; ET=ET_target)
r_target = view(x_target, 1:3)

every_day = 24 * 60 * 60

# NOTE: This is super brittle to initial conditions
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
# v⃗₀, v⃗ = lambert(; r⃗₀=x₀[1:3], r⃗=r_target, Δt, tₘ=1)
# x₀⁺ = vcat(x₀[1:3], v⃗₀)
# xₜ = vcat(r_target,v⃗) 
x₀⁺, xₜ = fixed_time_single_shoot(x₀, Δt, r_target; print_iter=false)

# Try coasting for a bit before starting terminal rendezvous
# x0 = [x_Earth[1:3];  -30404.024893932074; 8983.82902213017; 4411.591204327652]
x0 = x₀⁺
# x₀⁺ = x0

Δt_coast = 0.0*365*24*3600
x_coast = propagate_universal(x0, Δt_coast)

# hit an asteroid 
t0 = ET_start + Δt_coast
Δt_rendezvous = Δt - Δt_coast
# x_spacecraft, T_spacecraft, time_ET, mass_out_1 = calculate_rendezvous(x_coast, x_target, Δt_rendezvous, t0; m0=mining_ship.mass_total,
#                                                                     μ=GTOC12.μ_☉, dt=24*3600, output_times=every_day,
#                                                                     abstol = 1e-7, reltol=1e-9)

x_spacecraft, T_spacecraft, time_ET, mass_out_1, saved_values = calculate_rendezvous_from_Earth(x₀⁺, x_target, Δt, t0; m0=mining_ship.mass_total,
                                                                    μ=GTOC12.μ_☉, dt=3600, output_times=every_day,
                                                                    abstol = 1e-9, reltol=1e-11)

line_array = Vector{String}()
line_array, mining_ship = GTOC12.record_line!(line_array, "launch", [x_Earth, x_spacecraft[1]], ET_start, mining_ship)
                                                                    
ΔV_departure_1 = (x₀⁺-x_Earth)[4:6]
println("ΔV_departure_1 $(norm(ΔV_departure_1))")
r_burn_arc_1 = getindex.(x_spacecraft', 1:3)'
println("asteroid_rendezvous_resid $(x_spacecraft[end] - x_target)")


# TODO: shorten rendezvous arc to just be terminal:
# Therefore have less time while thrusting/disagreeing with their integration

println("mass pre burn 1: $(mining_ship.mass_total)")
line_array, mining_ship = GTOC12.record_line!(line_array, "burn", x_spacecraft, time_ET, mining_ship, control=saved_values.saveval, mass_burned=mass_out_1[1]-mass_out_1[end])
println("mass post burn 1: $(mining_ship.mass_total)")


#---------------------------------- Try plotting with our control too see if it makes sense -------------
function propagate_with_control!(ẋ, x, p, t)
    # unpack parameters
    interp_thrust_x, interp_thrust_y, interp_thrust_z, interp_mass = p
    thrust = [interp_thrust_x(t); interp_thrust_y(t); interp_thrust_z(t)]
    mass = interp_mass(t)

    # unpack state
    r⃗ = view(x, 1:3)
    v⃗ = view(x, 4:6)

    # Spacecraft Dynamics 
    # *************************************************
    a_grav = -GTOC12.μ_☉ * r⃗ / norm(r⃗)^3
    a_control = thrust/mass

    v̇ = a_grav + a_control
    ẋ[:] = [ v⃗; v̇ ]
    nothing
end

times = time_ET .- time_ET[1]
interp_thrust_x = LinearInterpolation(times, [val[1] for val in saved_values.saveval])
interp_thrust_y = LinearInterpolation(times, [val[2] for val in saved_values.saveval])
interp_thrust_z = LinearInterpolation(times, [val[3] for val in saved_values.saveval])
# interp_thrust_x = LinearInterpolation(times, [val for val in eachcol(T_spacecraft)[1]])
# interp_thrust_y = LinearInterpolation(times, [val for val in eachcol(T_spacecraft)[2]])
# interp_thrust_z = LinearInterpolation(times, [val for val in eachcol(T_spacecraft)[3]])
interp_mass   = LinearInterpolation(times, mass_out_1)
p = (interp_thrust_x, interp_thrust_y, interp_thrust_z, interp_mass)

# Pack up parameters and solve
prob = ODEProblem(propagate_with_control!, x_spacecraft[1], (times[1], times[end]), p)
sol2 = solve(prob, Vern9(), abstol = 1e-9, reltol=1e-11, tstops=every_day) 
x_sol = [state[1] for state in sol2.u]
y_sol = [state[2] for state in sol2.u]
z_sol = [state[3] for state in sol2.u]

# x_save_val = [state[1] for state in saved_values.saveval]
# y_save_val = [state[2] for state in saved_values.saveval]
# z_save_val = [state[3] for state in saved_values.saveval]
# plot(x_save_val)
# plot!(y_save_val)
# plot!(z_save_val)

# plot!(T_spacecraft, label="T_spacecraft")

#--------------------------------------------------------------------------------------------------------

# Deploy a miner
# TODO need to update mass in the mining_ship so that it can be used in this file 
# just output mining structure as well 
# could change this to event_line? 
# TODO mass is wrong in the second line when deploying 
# 
# hit an asteroid, deploy a miner
line_array, mining_ship = GTOC12.record_line!(line_array, "rendezvous", x_spacecraft[end], time_ET[end], mining_ship, rendez_flag="deploy", event_ID=ID_min)
println("mass post deploy 1: $(mining_ship.mass_total)")

# # Wait 10 years, then recover a miner
# t_wait = 10.5*365*24*3600
# ET_recover = 1.4776127999999104e9; #time_ET[end] + t_wait
# x_recover = propagate_universal(x_spacecraft[end], t_wait)
# line_array, mining_ship = GTOC12.record_line!(line_array, "rendezvous", x_recover, ET_recover, mining_ship, rendez_flag="recover", event_ID=ID_min)
# println("mass post recover 1: $(mining_ship.mass_total)")


# # Intercept the earth 
# ET_departure2 = ET_recover#0.25*365*24*3600
# transfer_time = 0.4*365*24*3600
# ET_rendevous2 = ET_departure2 + transfer_time
# transfer_time = ET_rendevous2 - ET_departure2
# x_intercept = get_body_state(Earth; ET=ET_rendevous2)
# x_burn_arc_2, T_vector, time_vector_ET, mass_out_2 = calculate_intercept(x_recover, x_intercept, transfer_time, ET_departure2; m0=mining_ship.mass_total, dt=24*3600, output_times=every_day,
# abstol=1e-9, 
# reltol=1e-11)
# r_burn_arc_2 = getindex.(x_burn_arc_2', 1:3)'

# println("earth_intercept_resid $(x_burn_arc_2[end] - x_intercept)")
# line_array, mining_ship = GTOC12.record_line!(line_array, "burn", x_burn_arc_2, time_vector_ET, mining_ship, control=T_vector, mass_burned=mass_out_2[1]-mass_out_2[end])
# println("mass post burn 2: $(mining_ship.mass_total)")

# # Record the Earth flyby and drop off cargo
# v_in_flyby = x_burn_arc_2[end][4:6]
# flyby_states = [vcat(x_intercept[1:3], v_in_flyby)]
# ET_flyby = time_vector_ET[end]
# # Calculate post-flyby state and append it to flyby_states vector
# x_flyby_out = naive_flyby(; x⃗_inrt=flyby_states[1], epoch_et=ET_flyby, flyby_body=Earth)
# push!(flyby_states, x_flyby_out)

# line_array, mining_ship = GTOC12.record_line!(line_array, "earth_flyby", flyby_states, ET_flyby, mining_ship)
# println("mass post Earth intercept: $(mining_ship.mass_total)")


# write lines to txt file
println("Printing results to file.")
file = open("Result.txt", "w")
for line in line_array
    write(file, line)
end
close(file)

# # Move file into problem directory
if isfile("problem/GTOC12_Verification_Program/GTOC12_Verification/Linux/Result.txt")
    rm("problem/GTOC12_Verification_Program/GTOC12_Verification/Linux/Result.txt")
end
mv("./Result.txt", "problem/GTOC12_Verification_Program/GTOC12_Verification/Linux/Result.txt")

# Run verification program
# run(``)

# Plotting
# Plot Earth/asteroid/spacecraft
ETs = [ET_start, ET_target]
plot_body(asteroid; ETs=ETs, color="Red")
plot_body!(GTOC12.Earth; ETs=ETs, color="Blue")
plot_coast!(x₀⁺, Δt; label="Coast 1", color="Green")
plot!(r_burn_arc_1[:,1], r_burn_arc_1[:,2], r_burn_arc_1[:,3], color="cyan", label="steered burn 1")
# plot!(r_burn_arc_2[:,1], r_burn_arc_2[:,2], r_burn_arc_2[:,3], color="magenta", label="steered burn 2")
plot!(x_sol, y_sol, z_sol, label = "extrapolated from controls", color="magenta")
