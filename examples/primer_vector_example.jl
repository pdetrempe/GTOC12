using GTOC12
using Plots
using LinearAlgebra
using ForwardDiff

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1))
asteroid = Asteroid(ID_min)

# ET of sim start (provide some margin for moving launch earlier)
ET_start = GTOC12.ET₀ + 365*24*3600

# Try out shooting method to hit asteroid
furnish_all_kernels()
DV₀ = [0; 0; 0]
x₀ = get_body_state(GTOC12.Earth; ET=ET_start) + [0; 0; 0; DV₀[:]]

# Transfer time
Δt = 1/4 * 365 * 24 * 3600
ΔV∞_max = 6000.0;

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = ET_start + Δt
x_target = get_body_state(asteroid; ET=ET_target)
r_target = view(x_target, 1:3)

# NOTE: This is super brittle to initial conditions
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
x₀⁺, xₜ = fixed_time_single_shoot(x₀, Δt, r_target; print_iter=false)
ΔV₀ = (x₀⁺-x₀)[4:6]
ΔVₜ = (xₜ-x_target)[4:6]

#--------------------- Plot Earth/asteroid/spacecraft
ETs = [ET_start, ET_target]
plot_body(asteroid; ETs=ETs, color=colormap("Reds"), ratio = 1 )
plot_body!(GTOC12.Earth; ETs=ETs, color=colormap("Blues"))
plot_coast!(x₀⁺, Δt; label="Coast 1", color=colormap("Greens"))

# Try to optimize burn by adding start/end coast
ΔV₀⁺, ΔVₜ⁺, t0, tf = optimize_impulsive_transfer(Earth, asteroid; ET_start=ET_start, Δt_guess=Δt, tol=1e-8)

# Plot optimized transfer
x_adjusted = get_body_state(GTOC12.Earth; ET=ET_start + t0) + [0; 0; 0; ΔV₀⁺[:]]
plot_coast!(x_adjusted, tf-t0; label="Adjusted Coast 1", color=colormap("Oranges"))

# # Plot primer vector as a function of transfer time
# times = range(Δt, 1.8*Δt, 100)
# p_vec = hcat([get_lambert_ΔV_and_p⃗(0, time; body_from=Earth, body_to=asteroid, ET_start=ET_start)[end] for time in times])


# Things to try next:
# - Non-dimensionalize/use canonical units
# - Try a joint Jacobian/Hessian (might need to use the cost function)
# - Just try a simple line search
