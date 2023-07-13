using DifferentialEquations

export calculate_rendezvous, calculate_intercept, low_thrust_optimal_control!, bc_rendezvous!, calculate_hamiltonian


function calculate_hamiltonian(x, p)
    m, λ, δ_star, u_star, T, μ = p
    # Solve entire problem in MEE
    A = A_equinoctial(x; μ=μ)
    B = B_equinoctial(x; μ=μ)

    H = T / m * δ_star + λ' * A + λ' * B * (T / m * δ_star * u_star)
    return H
end

function low_thrust_optimal_control!(dstate, state, p, t)
    # unpack parameters
    _, _, _, μ, CDU, CTU = p
    T = GTOC12.T_max/(CDU/CTU^2) # Max Thrust (kg-m/s^2). Non-dimensionalize

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
    u_star = -B' * λ / norm(B' * λ)
    S = norm(B' * λ) .- 1.0
    if S > 0
        δ_star = 1.0
    else
        δ_star = 0.0
    end
    cntrl = u_star*δ_star
    #println("Control: ", cntrl)

    # Spacecraft Dynamics 
    # *************************************************
    dx = A + T * δ_star / m * B * u_star

    # Mass variation
    # *************************************************
    # exhaust velocity
    g0 = GTOC12.g0/(CDU/CTU^2)    # m/s^2, non-dimensionalize
    Isp = canonical_time(GTOC12.Isp; CTU=CTU)  # seconds
    c = Isp * g0 # specific impulse and grav accel at sea level (m/s)
    dm = -T / c * δ_star

    # Costate diff eqs
    param = m, λ, δ_star, u_star, T, μ
    dH_dx = ForwardDiff.gradient(x -> calculate_hamiltonian(x, param), x)
    dλ = -dH_dx

    dstate[:] = [dx; dm; dλ]

end

function bc_rendezvous!(residual, state, p, t)
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

function bc_intercept!(residual, state, p, t)
    # u[1] is the beginning of the time span, and u[end] is the ending
    x0, m0, rf, μ, _, _ = p
    MEE_init = Cartesian2MEE(x0; μ=μ)
    MEE_current_canon = state
    MEE_current_end = state[end]
    x_end_canon = MEE2Cartesian(MEE_current_end[1:6]; μ=μ)

    # Initial/final boundary conditions
    p_0 = MEE_init
    p_f = rf

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
    residual[7] = x_end_canon[1] - p_f[1]
    residual[8] = x_end_canon[2] - p_f[2]
    residual[9] = x_end_canon[3] - p_f[3]

end


function calculate_rendezvous(x0, xf, Δt; m0, μ=GTOC12.μ_☉, dt=24*3600, abstol=1e-6, reltol=1e-10)
    # Largely a wrapper for the DifferentialEquations.jl 2-point BVP
    solve_bvp(x0, xf, Δt; boundary_condition=bc_rendezvous!, m0=m0, μ=μ, dt=dt, abstol=abstol, reltol=reltol)

end

function calculate_intercept(x0, xf, Δt; m0, μ=GTOC12.μ_☉, dt=24*3600, abstol=1e-6, reltol=1e-10)
    # Largely a wrapper for the DifferentialEquations.jl 2-point BVP
    solve_bvp(x0, xf, Δt; boundary_condition=bc_intercept!, m0=m0, μ=μ, dt=dt, abstol=abstol, reltol=reltol)

end

# Wrapper for boundary value problem with different conditions
function solve_bvp(x0, xf, Δt; boundary_condition, m0, μ=GTOC12.μ_☉, dt=24*3600, abstol=1e-6, reltol=1e-10)

    # Non-dimensionalize problem
    x0_canon, CDU, CTU, μ_canonical = get_canonical_state(x0; μ=μ)
    xf_canon = get_canonical_state(xf, CDU, CTU)
    MEE_init = Cartesian2MEE(x0_canon; μ=μ_canonical)
    
    state_init = vcat(MEE_init, m0, ones(6))
    tspan = canonical_time.((0.0, Δt); CTU=CTU)
    dt = canonical_time(dt; CTU=CTU)
    
    # Pack up parameters and solve
    p = (x0_canon, m0, xf_canon, μ_canonical, CDU, CTU)
    bvp2 = TwoPointBVProblem(low_thrust_optimal_control!, boundary_condition, state_init, tspan, p)
    sol = solve(bvp2, Shooting(Vern7()), dt=dt, abstol=abstol, reltol=reltol) # we need to use the MIRK4 solver for TwoPointBVProblem
    
    # Redimensionalize problem
    MEE_out = hcat([state[1:6] for state in sol.u])
    x_spacecraft = MEE2Cartesian.(MEE_out; μ=μ_canonical)
    redimensionalize_state!.(x_spacecraft; CDU=CDU, CTU=CTU)
    return x_spacecraft, sol
end




function calculation_post_process_burns(x0, sol)
    T_vector = zeros(length(sol.t), 3)
    time_vector_ET = zeros(length(sol.t))
    # Unpack Params from Non-dimensionalize problem
    μ = GTOC12.μ_☉
    #x0_canon, CDU, CTU, μ_canonical = get_canonical_state(sol.u[1][1:6]; μ=μ)
    for (idx, state_aug) in enumerate(sol.u)
        x = state_aug[1:6]
        CDU = norm(x[1:3])  # Canonical Distance Unit
        CTU = √(CDU^3/μ)    # Canonical Time Unit
        m = state_aug[7]
        λ = state_aug[8:13]
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
        T_vector[idx, :] = δ_star * u_star
        time_vector_ET[idx] = redimensionalize_time(sol.t[idx], CTU=CTU)
        # TODO update Mining Ship mass here?
    end
    # TODO convert T_vector to necessary units and frame 
    #time_vector_ET = canonical_time(sol.t, CTU=CTU)
    return T_vector, time_vector_ET
end

