using DifferentialEquations: TwoPointBVProblem, Shooting, Vern7, solve, SavingCallback, SavedValues

export calculate_rendezvous,
    calculate_rendezvous_from_Earth,
    calculate_intercept,
    low_thrust_optimal_control!,
    bc_rendezvous!,
    calculate_hamiltonian,
    calculation_post_process_burns

function calculate_hamiltonian(x, p)
    m, λ, δ_star, u_star, T, μ = p
    # Solve entire problem in MEE
    A = A_equinoctial(x; μ=μ)
    B = B_equinoctial(x; μ=μ)

    H = T / m * δ_star + λ' * A + λ' * B * (T / m * δ_star * u_star)
    return H
end

function save_control_rtn(state, t, integrator)
    _, _, _, μ, CDU, CTU = integrator.p

    # unpack state
    MEE = view(state, 1:6)
    λ = state[8:13]
    λp = λ[1]
    λp = redimensionalize_time(λp; CTU=1/CTU) # Redimensionalize only co-state with a length dimension (all others unitless)
    λ[1] = λp

    λv = view(λ, 2:6)
    redimensionalize_vel!(λv; CDU=CDU, CTU=CTU)


    x_cart = MEE2Cartesian(MEE; μ=μ)
    redimensionalize_state!(x_cart; CDU=CDU, CTU=CTU)
    MEE_dim = Cartesian2MEE(x_cart; μ=GTOC12.μ_☉)

    # Define B matrix (time varying)
    # *************************************************
    B = B_equinoctial(MEE_dim; μ=GTOC12.μ_☉)

    # Optimal Control Strategy 
    # *************************************************
    u_star = -B' * λ / norm(B' * λ)
    S = norm(B' * λ) .- 1.0
    if S > 0
        δ_star = 1.0
    else
        δ_star = 0.0
    end

    R_inrt2rtn = DCM_inertial_to_rtn(x_cart)

    return R_inrt2rtn * (GTOC12.T_max-10*eps()) * u_star * δ_star
end

function low_thrust_optimal_control!(dstate, state, p, t)
    # unpack parameters
    _, _, _, μ, CDU, CTU = p
    T = GTOC12.T_max / (CDU / CTU^2) # Max Thrust (kg-m/s^2). Non-dimensionalize

    # unpack state
    x = view(state, 1:6)
    m = state[7]
    λ = view(state, 8:13)

    # Define A and B matrices (time varying)
    # *************************************************
    # Solve entire problem in MEE
    A = A_equinoctial(x; μ=μ)
    B = B_equinoctial(x; μ=μ)

    # Optimal Control Strategy 
    # *************************************************
    u_star = -B' * λ / norm(B' * λ)
    S = norm(B' * λ) .- 1.0
    if S > 0
        δ_star = 1.0
    else
        δ_star = 0.0
    end
    cntrl = u_star * δ_star
    #println("Control: ", cntrl)

    # Spacecraft Dynamics 
    # *************************************************
    dx = A + T * δ_star / m * B * u_star

    # Mass variation
    # *************************************************
    # exhaust velocity
    g0 = GTOC12.g0 / (CDU / CTU^2)    # m/s^2, non-dimensionalize
    Isp = canonical_time(GTOC12.Isp; CTU=CTU)  # seconds
    c = Isp * g0 # specific impulse and grav accel at sea level (m/s)

    # Do ṁ calculations using dimensional units
    c_dim = GTOC12.Isp * GTOC12.g0
    dm = -T / c * δ_star


    # Costate diff eqs
    param = m, λ, δ_star, u_star, T, μ

    # Try speeding up using in-place ForwardDiff
    dH_dx = zeros(6) # where we'll be storing our results
    f = x -> calculate_hamiltonian(x, param)
    config = ForwardDiff.GradientConfig(f, x)
    ForwardDiff.gradient!(dH_dx, f, x, config)
    # dH_dx = ForwardDiff.gradient(x -> calculate_hamiltonian(x, param), x)
    dλ = -dH_dx

    dstate[1:6] = dx
    dstate[7] = dm
    dstate[8:13] = dλ
    nothing
end

function bc_rendezvous!(residual, state, p, t)
    # u[1] is the beginning of the time span, and u[end] is the ending
    x0, m0, xf, μ, CDU, CTU = p
    MEE_init = Cartesian2MEE(x0; μ=μ)
    MEE_target = Cartesian2MEE(xf; μ=μ)
    MEE_current_canon = state
    p_0 = MEE_init
    p_f = MEE_target

    pos2vel_ratio = 1e-4

    # Try converting to Cartesian for better approximation of actual Constraints
    MEE_current_end = state[end]
    x_end_canon = MEE2Cartesian(MEE_current_end[1:6]; μ=μ)
    x_start_canon = MEE2Cartesian(MEE_init; μ=μ)


    # initial boundary value 
    # ***********************************
    residual[1] = MEE_current_canon[1][1] - p_0[1]
    residual[2] = MEE_current_canon[1][2] - p_0[2]
    residual[3] = MEE_current_canon[1][3] - p_0[3]
    residual[4] = MEE_current_canon[1][4] - p_0[4]
    residual[5] = MEE_current_canon[1][5] - p_0[5]
    residual[6] = MEE_current_canon[1][6] - p_0[6]

    # final boundary value 
    # ***********************************
    residual[7] = MEE_current_canon[end][1] - p_f[1]
    residual[8] = MEE_current_canon[end][2] - p_f[2]
    residual[9] = MEE_current_canon[end][3] - p_f[3]
    residual[10] = MEE_current_canon[end][4] - p_f[4]
    residual[11] = MEE_current_canon[end][5] - p_f[5]
    residual[12] = MEE_current_canon[end][6] - p_f[6]

    # Initial mass constraint
    # ***********************************
    residual[13] = (state[1][7] - m0)

end

function bc_rendezvous_cartesian!(residual, state, p, t)
    # u[1] is the beginning of the time span, and u[end] is the ending
    x0, m0, xf, μ, CDU, CTU = p

    # Try converting to Cartesian for better approximation of actual Constraints
    MEE_current_end = state[end]
    MEE_current_start = state[1]
    x_end_canon = MEE2Cartesian(MEE_current_end[1:6]; μ=μ)
    x_start_canon = MEE2Cartesian(MEE_current_start[1:6]; μ=μ)

    pos_vel_ratio = 1e-6 # Constraint is 1000km vs. 1 m/s

    # initial boundary value 
    # ***********************************
    residual[1] = (x_start_canon[1] - x0[1]) * pos_vel_ratio * CDU
    residual[2] = (x_start_canon[2] - x0[2]) * pos_vel_ratio * CDU
    residual[3] = (x_start_canon[3] - x0[3]) * pos_vel_ratio * CDU
    residual[4] = (x_start_canon[4] - x0[4])*CDU/CTU
    residual[5] = (x_start_canon[5] - x0[5])*CDU/CTU
    residual[6] = (x_start_canon[6] - x0[6])*CDU/CTU

    # final boundary value 
    # ***********************************
    residual[7] = (x_end_canon[1] - xf[1]) * pos_vel_ratio * CDU
    residual[8] = (x_end_canon[2] - xf[2]) * pos_vel_ratio * CDU
    residual[9] = (x_end_canon[3] - xf[3]) * pos_vel_ratio * CDU
    residual[10] = (x_end_canon[4] - xf[4]) * CDU / CTU
    residual[11] = (x_end_canon[5] - xf[5]) * CDU / CTU
    residual[12] = (x_end_canon[6] - xf[6]) * CDU / CTU

    # Initial mass constraint
    # ***********************************
    residual[13] = (state[1][7] - m0)

end


function bc_rendezvous_cartesian_from_Earth!(residual, state, p, t)
    # For a rendezvous from Earth, we can actually vary the starting velocity and let the LV az/el clean up the rest
    x0, m0, xf, μ, CDU, CTU = p

    # Try converting to Cartesian for better approximation of actual Constraints
    MEE_current_end = state[end]
    MEE_current_start = state[1]
    x_end_canon = MEE2Cartesian(MEE_current_end[1:6]; μ=μ)
    x_start_canon = MEE2Cartesian(MEE_current_start[1:6]; μ=μ)

    pos_vel_ratio = 1e-6 # Constraint is 1000km vs. 1 m/s

    # initial boundary value 
    # ***********************************
    residual[1] = (x_start_canon[1] - x0[1]) * pos_vel_ratio * CDU
    residual[2] = (x_start_canon[2] - x0[2]) * pos_vel_ratio * CDU
    residual[3] = (x_start_canon[3] - x0[3]) * pos_vel_ratio * CDU
    # residual[4] = (x_start_canon[4] - x0[4])*CDU/CTU
    # residual[5] = (x_start_canon[5] - x0[5])*CDU/CTU
    # residual[6] = (x_start_canon[6] - x0[6])*CDU/CTU

    # final boundary value 
    # ***********************************
    residual[4] = (x_end_canon[1] - xf[1]) * pos_vel_ratio * CDU
    residual[5] = (x_end_canon[2] - xf[2]) * pos_vel_ratio * CDU
    residual[6] = (x_end_canon[3] - xf[3]) * pos_vel_ratio * CDU
    residual[7] = (x_end_canon[4] - xf[4]) * CDU / CTU
    residual[8] = (x_end_canon[5] - xf[5]) * CDU / CTU
    residual[9] = (x_end_canon[6] - xf[6]) * CDU / CTU

    # Initial mass constraint
    # ***********************************
    residual[10] = (state[1][7] - m0)

end

# TODO: Make Cartesian version of boundary condition
function bc_intercept_cartesian!(residual, state, p, t)
    x0, m0, xf, μ, CDU, CTU = p

    # Try converting to Cartesian for better approximation of actual Constraints
    MEE_current_end = state[end]
    MEE_current_start = state[1]
    x_end_canon = MEE2Cartesian(MEE_current_end[1:6]; μ=μ)
    x_start_canon = MEE2Cartesian(MEE_current_start[1:6]; μ=μ)

    pos_vel_ratio = 1e-6 # Constraint is 1000km vs. 1 m/s

    # initial boundary value 
    # ***********************************
    residual[1] = (x_start_canon[1] - x0[1]* CDU * pos_vel_ratio)
    residual[2] = (x_start_canon[2] - x0[2]* CDU * pos_vel_ratio)
    residual[3] = (x_start_canon[3] - x0[3]* CDU * pos_vel_ratio)
    residual[4] = (x_start_canon[4] - x0[4])* CDU/CTU
    residual[5] = (x_start_canon[5] - x0[5])* CDU/CTU
    residual[6] = (x_start_canon[6] - x0[6])* CDU/CTU

    # final boundary value 
    # ***********************************
    residual[7] = (x_end_canon[1] - xf[1]) * CDU * pos_vel_ratio
    residual[8] = (x_end_canon[2] - xf[2]) * CDU * pos_vel_ratio
    residual[9] = (x_end_canon[3] - xf[3]) * CDU * pos_vel_ratio

    # Initial mass constraint
    # ***********************************
    residual[10] = state[1][7] - m0
    nothing
end

function bc_intercept!(residual, state, p, t)
    # u[1] is the beginning of the time span, and u[end] is the ending
    x0, m0, rf, μ, _, _ = p
    MEE_init = Cartesian2MEE(x0; μ=μ)
    MEE_current_canon = state
    MEE_current_end = state[end]
    x_end_canon = MEE2Cartesian(MEE_current_end[1:6]; μ=μ)

    # Cartesian:MEE constraint ratio
    # Constraints are in 2 different state descriptions, need to tune to get similar performance
    cart2MEE_ratio = 10
    pos2vel_ratio = 1e-4

    # Initial/final boundary conditions
    p_0 = MEE_init
    p_f = rf

    # initial boundary value 
    # ***********************************
    residual[1] = (MEE_current_canon[1][1] - p_0[1])
    residual[2] = (MEE_current_canon[1][2] - p_0[2])
    residual[3] = (MEE_current_canon[1][3] - p_0[3])
    residual[4] = (MEE_current_canon[1][4] - p_0[4])
    residual[5] = (MEE_current_canon[1][5] - p_0[5])
    residual[6] = (MEE_current_canon[1][6] - p_0[6])

    # final boundary value 
    # ***********************************
    residual[7] = (x_end_canon[1] - p_f[1]) * cart2MEE_ratio
    residual[8] = (x_end_canon[2] - p_f[2]) * cart2MEE_ratio
    residual[9] = (x_end_canon[3] - p_f[3]) * cart2MEE_ratio

    # Initial mass constraint
    # ***********************************
    residual[13] = state[1][7] - m0
    nothing
end


function calculate_rendezvous(x0, xf, Δt, t0; m0, μ=GTOC12.μ_☉, dt=24 * 3600, abstol=1e-6, reltol=1e-10, output_times=nothing, kwargs...)
    # Largely a wrapper for the DifferentialEquations.jl 2-point BVP
    solve_bvp(x0, xf, Δt, t0; boundary_condition=bc_rendezvous_cartesian!, m0=m0, μ=μ, dt=dt, abstol=abstol, reltol=reltol, output_times, kwargs...)

end

function calculate_rendezvous_from_Earth(x0, xf, Δt, t0; m0, μ=GTOC12.μ_☉, dt=24 * 3600, abstol=1e-6, reltol=1e-10, output_times=nothing, kwargs...)
    # Largely a wrapper for the DifferentialEquations.jl 2-point BVP
    solve_bvp(x0, xf, Δt, t0; boundary_condition=bc_rendezvous_cartesian_from_Earth!, m0=m0, μ=μ, dt=dt, abstol=abstol, reltol=reltol, output_times, kwargs...)

end

function calculate_intercept(x0, xf, Δt, t0; m0, μ=GTOC12.μ_☉, dt=24 * 3600, abstol=1e-6, reltol=1e-10, output_times=nothing, kwargs...)
    # Largely a wrapper for the DifferentialEquations.jl 2-point BVP
    solve_bvp(x0, xf, Δt, t0; boundary_condition=bc_intercept_cartesian!, m0=m0, μ=μ, dt=dt, abstol=abstol, reltol=reltol, output_times, kwargs...)

end

# Wrapper for boundary value problem with different conditions
function solve_bvp(x0, xf, Δt, t0; boundary_condition, m0, μ=GTOC12.μ_☉, dt=24 * 3600, abstol=1e-6, reltol=1e-10, output_times=nothing, kwargs...)

    # Non-dimensionalize problem
    x0_canon, CDU, CTU, μ_canonical = get_canonical_state(x0; μ=μ)

    xf_canon = get_canonical_state(xf, CDU, CTU)
    MEE_init = Cartesian2MEE(x0_canon; μ=μ_canonical)
    state_init = vcat(MEE_init, m0, ones(6))
    tspan = canonical_time.((0.0, Δt); CTU=CTU)
    dt = canonical_time(dt; CTU=CTU)

    # Create a callback to save the control values as we go
    saved_values = SavedValues(Float64, Vector{Float64})

    # Pack up parameters and solve
    p = (x0_canon, m0, xf_canon, μ_canonical, CDU, CTU)
    bvp2 = TwoPointBVProblem(low_thrust_optimal_control!, boundary_condition, state_init, tspan, p)
    if output_times === nothing
        sol = solve(bvp2, Shooting(Vern9()), callback=cb, dt=dt, abstol=abstol, reltol=reltol, kwargs...) # we need to use the MIRK4 solver for TwoPointBVProblem
    else
        output_times = canonical_time(output_times; CTU=CTU)
        output_times_vec = [tspan[1]:output_times:tspan[2]; tspan[2]]
        cb = SavingCallback(save_control_rtn, saved_values; saveat=output_times_vec)
        sol = solve(bvp2, Shooting(Vern9()), callback=cb, dt=dt, abstol=abstol, reltol=reltol, saveat=output_times, kwargs...) # we need to use the MIRK4 solver for TwoPointBVProblem
    end

    T_vector, time_vector_ET = calculate_optimal_control(sol, t0, CDU=CDU, CTU=CTU, μ=μ_canonical)

    # Redimensionalize problem
    MEE_out = hcat([state[1:6] for state in sol.u])
    mass_out = [state[7] for state in sol.u]
    x_spacecraft = MEE2Cartesian.(MEE_out; μ=μ_canonical)
    redimensionalize_state!.(x_spacecraft; CDU=CDU, CTU=CTU)
    # saved_values.t *= CTU # NOTE: THis doesn't work, so beware of time on saved_values
    return x_spacecraft, T_vector, time_vector_ET, mass_out, saved_values
end

function calculate_optimal_control(sol, t0; CDU, CTU, μ)
    T = (GTOC12.T_max - 10*eps())/ (CDU / CTU^2) # Max Thrust (kg-m/s^2). Non-dimensionalize
    T_vector = zeros(length(sol.t), 3)
    time_vector_ET = zeros(length(sol.t))
    for (idx, state_aug) in enumerate(sol.u)
        MEE = state_aug[1:6]
        m = state_aug[7]
        λ = state_aug[8:13]
        B = B_equinoctial(MEE; μ=μ)

        # Optimal Control Strategy 
        # *************************************************
        u_star = -B' * λ / norm(B' * λ)
        S = norm(B' * λ) .- 1.0
        if S > 0
            δ_star = 1.0
        else
            δ_star = 0.0
        end
        # Need to rotate control from RTN frame to inertial frame
        x_cart = MEE2Cartesian(MEE; μ=μ)
        R_inrt2rtn = DCM_inertial_to_rtn(x_cart)
        T_vector[idx, :] = R_inrt2rtn' * T * δ_star * u_star
        time_vector_ET[idx] = redimensionalize_time(sol.t[idx], CTU=CTU)
        # TODO update Mining Ship mass here?
    end
    # Convert T_vector to necessary units
    T_vector *= (CDU / CTU^2)
    # add t0 to delta t vector
    time_vector_ET .+= t0
    #time_vector_ET = canonical_time(sol.t, CTU=CTU)
    return T_vector, time_vector_ET
end

