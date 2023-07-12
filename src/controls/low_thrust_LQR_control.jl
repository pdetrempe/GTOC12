using GTOC12
using Plots
using LinearAlgebra
using DifferentialEquations
using ForwardDiff
# using BoundaryValueDiffEq

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1.0))
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

m = 500.0 # kg

# Non-dimensionalize problem
x0_canon, CDU, CTU, μ_canonical = get_canonical_state(x₀⁺; μ=GTOC12.μ_☉)
xf_canon = get_canonical_state(x_target, CDU, CTU)
MEE_init = Cartesian2MEE(x0_canon; μ=μ_canonical)

state_init = vcat(MEE_init, m, ones(6))
tspan = canonical_time.((0.0, Δt); CTU=CTU)
dt = canonical_time(24 * 3600.0; CTU=CTU)

# Pack up parameters and solve
p = (x0_canon, m, xf_canon, μ_canonical, CDU, CTU)
bvp2 = TwoPointBVProblem(low_thrust_optimal_control!, bc2!, state_init, tspan, p)
sol = solve(bvp2, Shooting(Vern7()), dt=dt, abstol=1e-6, reltol=1e-10) # we need to use the MIRK4 solver for TwoPointBVProblem

# Redimensionalize problem
MEE_out = hcat([state[1:6] for state in sol.u])
x_spacecraft = MEE2Cartesian.(MEE_out; μ=μ_canonical)
redimensionalize_state!.(x_spacecraft; CDU=CDU, CTU=CTU)
r_spacecraft = getindex.(x_spacecraft', 1:3)'

# Plot solution alongside coasts
plot!(r_spacecraft[:,1], r_spacecraft[:,2], r_spacecraft[:,3], color="cyan", label="steered burn")

# TODO: Figure out how to use this converged solution to feed to the optimizer
