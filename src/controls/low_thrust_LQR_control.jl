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

function low_thrust_optimal_control!(dstate, state, p, t)
    # unpack parameters
    _, _, _, μ, CDU, CTU = p
    T = 0.6 /(CDU/CTU^2) # Max Thrust (kg-m/s^2). Non-dimensionalize

    # unpack state
    x = state[1:6]
    m = state[7]
    λ = state[8:13]

    # Define A and B matrices (time varying)
    # *************************************************
    # Solve entire problem in MEE
    A = A_equinoctial(x; μ=μ)
    B = B_equinoctial(x; μ=μ)

    # Optimal Control Strategy 
    # *************************************************
    #u_star = [0,0,0]
    #δ_star = 0
    u_star = -B' * λ / norm(B' * λ)
    S = norm(B' * λ) .- 1.0
    if S > 0
        δ_star = 1.0
    else
        δ_star = 0.0
    end

    # Spacecraft Dynamics 
    # *************************************************
    dx = A + T * δ_star / m * B * u_star

    # Mass variation
    # *************************************************
    # exhaust velocity
    g0 = 9.81/(CDU/CTU^2)    # m/s^2, non-dimensionalize
    Isp = canonical_time(4000.0; CTU=CTU)  # seconds
    c = Isp * g0 # specific impulse and grav accel at sea level (m/s)
    dm = -T / c * δ_star

    # Costate diff eqs
    param = m, λ, δ_star, u_star, T, μ
    dH_dx = ForwardDiff.gradient(x -> calculate_hamiltonian(x, param), x) #[1:6])
    dλ = -dH_dx

    dstate[:] = [dx; dm; dλ]

end

function calculate_hamiltonian(x, p)
    m, λ, δ_star, u_star, T, μ = p
    # Solve entire problem in MEE
    A = A_equinoctial(x; μ=μ)
    B = B_equinoctial(x; μ=μ)

    H = T / m * δ_star + λ' * A + λ' * B * (T / m * δ_star * u_star)
    return H
end

function bc2!(residual, state, p, t)
    # u[1] is the beginning of the time span, and u[end] is the ending
    x0, m0, xf, μ, _, _ = p
    MEE_init = Cartesian2MEE(x0; μ=μ)
    MEE_target = Cartesian2MEE(xf; μ=μ)
    MEE_current_canon = state
    p_0 = MEE_init
    p_f = MEE_target

    # initial boundary value 
    # ***********************************
    residual[1] = MEE_current_canon[1][1] - p_0[1]
    residual[2] = MEE_current_canon[1][2] - p_0[2]
    residual[3] = MEE_current_canon[1][3] - p_0[3]
    residual[4] = MEE_current_canon[1][4] - p_0[4]
    residual[5] = MEE_current_canon[1][5] - p_0[5]
    residual[6] = MEE_current_canon[1][6] - p_0[6]
    residual[7] = state[1][7] - m0 # This mass IC doesn't seem to affect the outcome

    # final boundary value 
    # ***********************************
    residual[7] = MEE_current_canon[end][1] - p_f[1]
    residual[8] = MEE_current_canon[end][2] - p_f[2]
    residual[9] = MEE_current_canon[end][3] - p_f[3]
    residual[10] = MEE_current_canon[end][4] - p_f[4]
    residual[11] = MEE_current_canon[end][5] - p_f[5]
    residual[12] = MEE_current_canon[end][6] - p_f[6]
end

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
