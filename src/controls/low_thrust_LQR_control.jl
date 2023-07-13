using GTOC12
using Plots
using LinearAlgebra


# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1.0))
# ID_min = 38
asteroid = Asteroid(ID_min)

# Try out shooting method to hit asteroid
furnish_all_kernels()
DV₀ = [-1000.0; 0.0; 3500.0]
x₀ = get_body_state(GTOC12.Earth; ET=GTOC12.ET₀) + [0; 0; 0; DV₀[:]]

# Transfer time
Δt = 0.75 * 365 * 24 * 3600
ΔV∞_max = 6000.0;

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = GTOC12.ET₀ + Δt
x_target = get_body_state(asteroid; ET=ET_target)
r_target = view(x_target, 1:3)

# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
x₀⁺, xₜ = fixed_time_single_shoot(x₀, Δt, r_target; print_iter=false)
ΔV_departure_1 = (x₀⁺-x₀)[4:6]
ΔV_arrival_1 = (xₜ-x_target)[4:6]

# Plot Earth/asteroid/spacecraft
ETs = [GTOC12.ET₀, ET_target]
plot_body(asteroid; ETs=ETs, color="Red")
plot_body!(GTOC12.Earth; ETs=ETs, color="Blue")
plot_coast!(x₀⁺, Δt; label="Coast 1", color="Green")

# TODO: Optimize the above for DV/via Time-of-flight

# TODO: For launch case, just need initial time, not states
#states_out, controls_out, time_out = optimize_impulsive_launch(x₀, x_target, Δt; ΔV₀=ΔV_departure_1, asteroid_ID=ID_min)

m = 3000.0 # kg

# TODO:
# - Extract burn directions (in inertial space)
# - Extract times (specify times to save at, e.g. daily)
# - Extract mass

@time x_burn_arc_1, sol = calculate_rendezvous(x₀⁺, x_target, Δt; m0=m, μ=GTOC12.μ_☉, dt=24*3600,
abstol=1e-5, 
reltol=1e-8)
r_burn_arc_1 = getindex.(x_burn_arc_1', 1:3)'
mass_arc_1 = getindex.(sol.u, 7)

# # # # Calculate re-rendezvous with asteroid years later
# ET_departure2 = ET_target + 0#0.25*365*24*3600
# transfer_time = 0.35*365*24*3600
# ET_rendevous2 = ET_departure2 + transfer_time
# transfer_time = ET_rendevous2 - ET_departure2
# x_rendezvous2 = get_body_state(Earth; ET=ET_rendevous2)
# x0 = get_body_state(asteroid; ET=ET_departure2)

# # The BVP will try to use the final value as an initial guess, which is bad when we just want to intercept
# v⃗₀_lambert, vf_lambert = lambert(; r⃗₀=x0[1:3], r⃗=x_rendezvous2[1:3], Δt=transfer_time, tₘ=-1)
# xf = [x_rendezvous2[1:3]; vf_lambert]

# # Try intercepting this Lambert arc partway around
# x0_lambert = [x0[1:3]; v⃗₀_lambert]
# t_to_lambert = 1*transfer_time
# x_lambert_some_time_around = propagate_universal(x0_lambert, t_to_lambert)
# x_burn_arc_2 = calculate_intercept(x0, x_lambert_some_time_around, t_to_lambert; m0=m, dt=24*3600)
# r_burn_arc_2 = getindex.(x_burn_arc_2', 1:3)'

# V_∞ = norm((x_burn_arc_2[end] - x_rendezvous2)[4:6])
# println("V_∞ = $V_∞")

# # Plot solution alongside coasts
# plot!(r_burn_arc_1[:,1], r_burn_arc_1[:,2], r_burn_arc_1[:,3], color="cyan", label="steered burn 1")
# plot!(r_burn_arc_2[:,1], r_burn_arc_2[:,2], r_burn_arc_2[:,3], color="magenta", label="steered burn 2")

# plot_coast!(x0_lambert, transfer_time; label="Lambert arc 2", color="Green")
# plot_body!(Earth; ETs=[ET_departure2, ET_rendevous2], color="Cyan")

# TODO: Figure out how to use this converged solution to feed to the optimizer
