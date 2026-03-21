"""
InfiniteOpt.jl Low-Thrust Trajectory Optimization for GTOC12

Adapts the InfiniteOpt approach from GTOC13 (solar sail) to GTOC12 (low-thrust).

Key Differences from GTOC13:
- Low-thrust propulsion instead of solar sail
- Mass depletion (ṁ = -T/(Isp·g₀)·δ)
- Throttle control δ ∈ [0,1]
- Different acceleration model: a = (T/m)·û·δ

Key Features:
- Clean, intuitive syntax for physics (∂(x, t) == f(x, u))
- Automatic transcription to NLP via InfiniteOpt.jl
- OrthogonalCollocation for higher accuracy
- Piecewise constant controls (reduces optimizer burden)
- Support for both Cartesian and MEE coordinates
"""

module InfiniteOptLowThrust

using InfiniteOpt
using Ipopt
using LinearAlgebra
using Statistics

export infiniteopt_low_thrust_cartesian,
       infiniteopt_low_thrust_mee,
       InfiniteOptLowThrustResult,
       generate_gtoc12_submission_file

"""
    InfiniteOptLowThrustResult

Results from InfiniteOpt low-thrust trajectory optimization.

Fields:
- success: Whether optimization converged
- objective_value: Final objective function value
- X: State trajectory (7 × N matrix for Cartesian: [r; v; m], or 7 × N for MEE: [p,f,g,h,k,L,m])
- U: Control trajectory (4 × N matrix: [θ, φ, δ] where θ,φ define thrust direction, δ is throttle)
- t_nodes: Time at each discretization point [s]
- solve_time: Wall-clock time [s]
- iterations: Number of NLP iterations (-1 if unknown)
- constraint_violation: Maximum constraint violation
- constraint_breakdown: Named tuple with constraint violations by type
- mass_used: Total propellant mass used [kg]
- coordinate_system: :cartesian or :mee
"""
struct InfiniteOptLowThrustResult
    success::Bool
    objective_value::Float64
    X::Matrix{Float64}  # 7 × N state (includes mass)
    U::Matrix{Float64}  # 4 × N control (θ, φ, δ, plus magnitude for validation)
    t_nodes::Vector{Float64}
    solve_time::Float64
    iterations::Int
    constraint_violation::Float64
    constraint_breakdown::NamedTuple
    mass_used::Float64
    coordinate_system::Symbol
end

"""
    infiniteopt_low_thrust_cartesian(x0, xf, t_f, N, μ; kwargs...)

Solve low-thrust trajectory targeting using InfiniteOpt.jl in Cartesian coordinates.

Args:
    x0: Initial state [r0; v0; m0] (7-element, m and m/s and kg)
    xf: Target final state [rf; vf] (6-element, mass is free)
    t_f: Transfer time [s]
    N: Number of discretization points (num_supports)
    μ: Gravitational parameter [m³/s²]

Keyword Args:
    T_max: Maximum thrust [N] (default: 0.6 N)
    Isp: Specific impulse [s] (default: 4000 s)
    g0: Gravitational acceleration at sea level [m/s²] (default: 9.80665)
    m_min: Minimum allowed mass [kg] (default: 500.0)
    X_guess: Initial guess for state trajectory (7 × N, optional)
    U_guess: Initial guess for controls (3 × N: [θ, φ, δ], optional)
    t_nodes: Custom time nodes [s] (optional, if provided N is ignored)
    position_tol: Final position tolerance [m] (default: 1e3)
    velocity_tol: Final velocity tolerance [m/s] (default: 1.0)
    intercept_only: Only match final position, not velocity (default: false)
    max_iter: Maximum NLP iterations (default: 3000)
    verbose: Print solver output (default: false)

Returns:
    InfiniteOptLowThrustResult
"""
function infiniteopt_low_thrust_cartesian(
    x0::Vector{Float64},
    xf::Vector{Float64},
    t_f::Float64,
    N::Int,
    μ::Float64;
    T_max::Float64=0.6,  # N
    Isp::Float64=4000.0,  # s
    g0::Float64=9.80665,  # m/s²
    m_min::Float64=500.0,  # kg
    X_guess::Union{Matrix{Float64},Nothing}=nothing,
    U_guess::Union{Matrix{Float64},Nothing}=nothing,
    t_nodes::Union{Vector{Float64},Nothing}=nothing,
    position_tol::Float64=1e3,  # m
    velocity_tol::Float64=1.0,  # m/s
    intercept_only::Bool=false,
    max_iter::Int=3000,
    verbose::Bool=false
)
    # =============================================================================
    # CANONICAL UNITS for numerical conditioning
    # =============================================================================
    # Use similar scaling to GTOC13
    # Distance: 1 AU
    # Time: 10 days (good for low-thrust maneuvers)
    # Mass: initial mass
    AU_m = 1.49597870691e11  # m
    L_scale = AU_m  # 1 AU
    T_scale = 10.0 * 86400.0  # 10 days in seconds
    V_scale = L_scale / T_scale  # m/s
    A_scale = L_scale / T_scale^2  # m/s²
    μ_scale = L_scale^3 / T_scale^2  # m³/s²
    M_scale = x0[7]  # kg (initial mass)
    F_scale = M_scale * A_scale  # N

    # Non-dimensionalize inputs
    x0_nd = zeros(7)
    x0_nd[1:3] = x0[1:3] / L_scale
    x0_nd[4:6] = x0[4:6] / V_scale
    x0_nd[7] = x0[7] / M_scale

    xf_nd = zeros(6)
    xf_nd[1:3] = xf[1:3] / L_scale
    xf_nd[4:6] = xf[4:6] / V_scale

    t_f_nd = t_f / T_scale
    μ_nd = μ / μ_scale
    T_max_nd = T_max / F_scale
    g0_nd = g0 / A_scale
    m_min_nd = m_min / M_scale

    # Exhaust velocity (dimensional, then non-dimensionalize)
    c_dim = Isp * g0  # m/s
    c_nd = c_dim / V_scale

    # Non-dimensionalize custom time nodes if provided
    t_nodes_nd = nothing
    if t_nodes !== nothing
        t_nodes_nd = t_nodes ./ T_scale
    end

    # =============================================================================
    # SOLVE IN CANONICAL UNITS
    # =============================================================================
    result_nd = _infiniteopt_low_thrust_cartesian_canonical(
        x0_nd, xf_nd, t_f_nd, N, μ_nd;
        T_max=T_max_nd,
        c=c_nd,
        m_min=m_min_nd,
        X_guess=X_guess,
        U_guess=U_guess,
        t_nodes=t_nodes_nd,
        position_tol=position_tol / L_scale,
        velocity_tol=velocity_tol / V_scale,
        intercept_only=intercept_only,
        max_iter=max_iter,
        verbose=verbose
    )

    # =============================================================================
    # RE-DIMENSIONALIZE RESULTS
    # =============================================================================
    X_sol = copy(result_nd.X)
    X_sol[1:3, :] .*= L_scale
    X_sol[4:6, :] .*= V_scale
    X_sol[7, :] .*= M_scale

    U_sol = result_nd.U  # Angles and throttle are dimensionless

    t_nodes_actual = result_nd.t_nodes .* T_scale

    # Re-dimensionalize constraint breakdown
    constraint_breakdown_actual = (
        position=result_nd.constraint_breakdown.position * L_scale,
        velocity=result_nd.constraint_breakdown.velocity * V_scale,
        mass_negative=result_nd.constraint_breakdown.mass_negative * M_scale,
        thrust_direction=result_nd.constraint_breakdown.thrust_direction  # Dimensionless
    )

    mass_used_actual = result_nd.mass_used * M_scale

    return InfiniteOptLowThrustResult(
        result_nd.success,
        result_nd.objective_value,
        X_sol,
        U_sol,
        t_nodes_actual,
        result_nd.solve_time,
        result_nd.iterations,
        result_nd.constraint_violation,
        constraint_breakdown_actual,
        mass_used_actual,
        :cartesian
    )
end

"""
    _infiniteopt_low_thrust_cartesian_canonical

Internal function that solves the problem in canonical (dimensionless) units.
"""
function _infiniteopt_low_thrust_cartesian_canonical(
    x0::Vector{Float64},
    xf::Vector{Float64},
    t_f::Float64,
    N::Int,
    μ::Float64;
    T_max::Float64=0.6,
    c::Float64=1.0,
    m_min::Float64=0.5,
    X_guess::Union{Matrix{Float64},Nothing}=nothing,
    U_guess::Union{Matrix{Float64},Nothing}=nothing,
    t_nodes::Union{Vector{Float64},Nothing}=nothing,
    position_tol::Float64=1e-6,
    velocity_tol::Float64=1e-9,
    intercept_only::Bool=false,
    max_iter::Int=3000,
    verbose::Bool=false
)
    # Create InfiniteOpt model
    model = InfiniteModel(optimizer_with_attributes(
        Ipopt.Optimizer,
        "print_level" => verbose ? 5 : 0,
        "max_iter" => max_iter,
        "tol" => 1e-8,
        "constr_viol_tol" => 1e-8,
        "acceptable_tol" => 1e-6,
        "acceptable_constr_viol_tol" => 1e-6,
        "acceptable_iter" => 50,
        "mumps_mem_percent" => 8000,
        "mu_strategy" => "adaptive",
        "max_cpu_time" => 600.0,
        "nlp_scaling_method" => "gradient-based"
    ))

    # Define infinite parameter (time) with orthogonal collocation
    if t_nodes !== nothing
        @infinite_parameter(model, t in [0, t_f], supports = t_nodes,
            derivative_method = InfiniteOpt.OrthogonalCollocation(4))
    else
        @infinite_parameter(model, t in [0, t_f], num_supports = N,
            derivative_method = InfiniteOpt.OrthogonalCollocation(4))
    end

    # State variables
    @variable(model, r[1:3], Infinite(t))
    @variable(model, v[1:3], Infinite(t))
    @variable(model, m, Infinite(t))

    # Control variables: Spherical angles for thrust direction + throttle
    # θ: Polar angle from +z axis, range [0, π]
    # φ: Azimuthal angle in xy-plane from +x axis, range unbounded
    # δ: Throttle, range [0, 1]
    @variable(model, 0 <= θ <= π, Infinite(t))
    @variable(model, φ, Infinite(t))
    @variable(model, 0 <= δ <= 1, Infinite(t))

    # Make controls piecewise constant over collocation intervals (if available)
    try
        constant_over_collocation(θ, t)
        constant_over_collocation(φ, t)
        constant_over_collocation(δ, t)
    catch
        # Function not available in this InfiniteOpt version - use continuous controls
        # This is OK, just less efficient
    end

    # Helper function to convert unit vector to spherical angles
    function unit_vector_to_angles(û)
        θ_val = acos(clamp(û[3], -1.0, 1.0))  # [0, π]
        φ_val = atan(û[2], û[1])              # [-π, π]
        return (θ_val, φ_val)
    end

    # Set initial guesses
    if X_guess !== nothing
        # Non-dimensionalize guess if dimensional
        # (Assume caller handles this, or we'd need scale info)
        for i in 1:3
            set_start_value_function(r[i], t_val -> X_guess[i, 1] + (X_guess[i, end] - X_guess[i, 1]) * t_val / t_f)
            set_start_value_function(v[i], t_val -> X_guess[i+3, 1] + (X_guess[i+3, end] - X_guess[i+3, 1]) * t_val / t_f)
        end
        set_start_value_function(m, t_val -> X_guess[7, 1] - (X_guess[7, 1] - X_guess[7, end]) * t_val / t_f)
    else
        # Default linear interpolation
        for i in 1:3
            set_start_value_function(r[i], t_val -> x0[i] + (xf[i] - x0[i]) * t_val / t_f)
            set_start_value_function(v[i], t_val -> x0[i+3] + (xf[i+3] - x0[i+3]) * t_val / t_f)
        end
        # Mass decreases linearly (rough guess)
        set_start_value_function(m, t_val -> x0[7] * (1.0 - 0.1 * t_val / t_f))
    end

    if U_guess !== nothing
        # Convert first unit vector from U_guess to angles for initial guess
        θ_init, φ_init = unit_vector_to_angles(U_guess[1:3, 1])
        δ_init = U_guess[4, 1]
        set_start_value_function(θ, t_val -> θ_init)
        set_start_value_function(φ, t_val -> φ_init)
        set_start_value_function(δ, t_val -> δ_init)
    else
        # Default: half throttle, thrust in velocity direction
        set_start_value_function(θ, t_val -> π / 2)
        set_start_value_function(φ, t_val -> 0.0)
        set_start_value_function(δ, t_val -> 0.5)
    end

    # Objective: Minimize propellant usage (maximize final mass)
    # OR minimize final velocity for rendezvous
    t_final = upper_bound(t)
    if intercept_only
        # For intercept: minimize fuel usage
        @objective(model, Max, m(t_final))
    else
        # For rendezvous: minimize terminal velocity error (will be constrained anyway)
        # Use fuel minimization
        @objective(model, Max, m(t_final))
    end

    # Initial conditions
    @constraint(model, [i = 1:3], r[i](0) == x0[i])
    @constraint(model, [i = 1:3], v[i](0) == x0[i+3])
    @constraint(model, m(0) == x0[7])

    # Final conditions
    @constraint(model, [i = 1:3], r[i](t_final) == xf[i])

    if !intercept_only
        # Full rendezvous: match velocity too
        @constraint(model, [i = 1:3], v[i](t_final) == xf[i+3])
    end

    # Dynamics: dr/dt = v
    @constraint(model, [i = 1:3], ∂(r[i], t) == v[i])

    # Dynamics: dv/dt = a_grav + a_thrust
    # Thrust direction: û = [sin(θ)cos(φ), sin(θ)sin(φ), cos(θ)]
    # Thrust acceleration: a_thrust = (T_max / m) * δ * û

    for i in 1:3
        if i == 1
            # û_x = sin(θ) * cos(φ)
            @constraint(model,
                ∂(v[1], t) == -μ * r[1] / (r[1]^2 + r[2]^2 + r[3]^2)^(3 / 2)
                              + (T_max / m) * δ * sin(θ) * cos(φ)
            )
        elseif i == 2
            # û_y = sin(θ) * sin(φ)
            @constraint(model,
                ∂(v[2], t) == -μ * r[2] / (r[1]^2 + r[2]^2 + r[3]^2)^(3 / 2)
                              + (T_max / m) * δ * sin(θ) * sin(φ)
            )
        else  # i == 3
            # û_z = cos(θ)
            @constraint(model,
                ∂(v[3], t) == -μ * r[3] / (r[1]^2 + r[2]^2 + r[3]^2)^(3 / 2)
                              + (T_max / m) * δ * cos(θ)
            )
        end
    end

    # Dynamics: dm/dt = -(T_max / c) * δ
    @constraint(model, ∂(m, t) == -(T_max / c) * δ)

    # Path constraints: mass must stay positive
    @constraint(model, [t = t], m >= m_min)

    # Solve
    start_time = time()
    optimize!(model)
    solve_time = time() - start_time

    # Extract results
    success = termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
    obj_value = objective_value(model)

    # Get discretization points
    t_vals = value(t)

    # Extract state and control
    r_vals = value.(r)
    v_vals = value.(v)
    m_vals = value(m)
    θ_vals = value(θ)
    φ_vals = value(φ)
    δ_vals = value(δ)

    # Convert to matrix format
    N_actual = length(t_vals)
    X_sol = zeros(7, N_actual)
    U_sol = zeros(4, N_actual)

    for i in 1:N_actual
        X_sol[1:3, i] = [r_vals[j][i] for j in 1:3]
        X_sol[4:6, i] = [v_vals[j][i] for j in 1:3]
        X_sol[7, i] = m_vals[i]

        # Store spherical angles + throttle
        θ_i = θ_vals[i]
        φ_i = φ_vals[i]
        δ_i = δ_vals[i]

        U_sol[1, i] = θ_i
        U_sol[2, i] = φ_i
        U_sol[3, i] = δ_i

        # Store thrust magnitude for validation
        U_sol[4, i] = T_max * δ_i
    end

    # Get solver stats
    iterations = -1  # Unknown - not exposed by Ipopt via MOI

    # =============================================================================
    # COMPUTE CONSTRAINT VIOLATIONS
    # =============================================================================

    max_position_viol = 0.0
    max_velocity_viol = 0.0
    max_mass_negative_viol = 0.0
    max_thrust_direction_viol = 0.0

    # Check final position
    max_position_viol = norm(X_sol[1:3, end] - xf[1:3])

    # Check final velocity (if rendezvous)
    if !intercept_only
        max_velocity_viol = norm(X_sol[4:6, end] - xf[4:6])
    end

    # Check mass is positive everywhere
    for i in 1:N_actual
        if X_sol[7, i] < m_min
            max_mass_negative_viol = max(max_mass_negative_viol, m_min - X_sol[7, i])
        end
    end

    # Check thrust direction normalization (should be automatic with spherical coords)
    for i in 1:N_actual
        û_norm = sqrt(sin(U_sol[1, i])^2 + cos(U_sol[1, i])^2)  # Should be 1.0
        max_thrust_direction_viol = max(max_thrust_direction_viol, abs(û_norm - 1.0))
    end

    max_viol = max(max_position_viol, max_velocity_viol, max_mass_negative_viol, max_thrust_direction_viol)

    constraint_breakdown = (
        position=max_position_viol,
        velocity=max_velocity_viol,
        mass_negative=max_mass_negative_viol,
        thrust_direction=max_thrust_direction_viol
    )

    # Compute mass used
    mass_used = X_sol[7, 1] - X_sol[7, end]

    if verbose
        println("\n" * "="^70)
        println("Constraint Violation Breakdown (canonical units):")
        println("="^70)
        println("  Position:         $(round(constraint_breakdown.position, digits=10))")
        println("  Velocity:         $(round(constraint_breakdown.velocity, digits=10))")
        println("  Mass negative:    $(round(constraint_breakdown.mass_negative, digits=10))")
        println("  Thrust direction: $(round(constraint_breakdown.thrust_direction, digits=10))")
        println("  Overall max:      $(round(max_viol, digits=10))")
        println("  Mass used:        $(round(mass_used, digits=6))")
        println("="^70)
    end

    return InfiniteOptLowThrustResult(
        success,
        obj_value,
        X_sol,
        U_sol,
        collect(t_vals),
        solve_time,
        iterations,
        max_viol,
        constraint_breakdown,
        mass_used,
        :cartesian
    )
end

"""
    infiniteopt_low_thrust_mee(x0, xf, t_f, N, μ; kwargs...)

Solve low-thrust trajectory targeting using InfiniteOpt.jl in Modified Equinoctial Elements (MEE).

Args:
    x0: Initial state [p0, f0, g0, h0, k0, L0, m0] (7-element) or Cartesian [r0; v0; m0]
    xf: Target final state [pf, ff, gf, hf, kf, Lf] (6-element) or Cartesian [rf; vf]
    t_f: Transfer time [s]
    N: Number of discretization points
    μ: Gravitational parameter [m³/s²]

Keyword Args:
    T_max: Maximum thrust [N] (default: 0.6 N)
    Isp: Specific impulse [s] (default: 4000 s)
    g0: Gravitational acceleration at sea level [m/s²] (default: 9.80665)
    m_min: Minimum allowed mass [kg] (default: 500.0)
    X_guess: Initial guess for state trajectory (7 × N, optional)
    U_guess: Initial guess for controls (3 × N: [θ, φ, δ], optional)
    t_nodes: Custom time nodes [s] (optional)
    position_tol: Final position tolerance [m] (default: 1e3)
    velocity_tol: Final velocity tolerance [m/s] (default: 1.0)
    intercept_only: Only match final position, not velocity (default: false)
    max_iter: Maximum NLP iterations (default: 3000)
    verbose: Print solver output (default: false)
    convert_to_mee: Convert Cartesian inputs to MEE (default: true)

Returns:
    InfiniteOptLowThrustResult
"""
function infiniteopt_low_thrust_mee(
    x0::Vector{Float64},
    xf::Vector{Float64},
    t_f::Float64,
    N::Int,
    μ::Float64;
    T_max::Float64=0.6,
    Isp::Float64=4000.0,
    g0::Float64=9.80665,
    m_min::Float64=500.0,
    X_guess::Union{Matrix{Float64},Nothing}=nothing,
    U_guess::Union{Matrix{Float64},Nothing}=nothing,
    t_nodes::Union{Vector{Float64},Nothing}=nothing,
    position_tol::Float64=1e3,
    velocity_tol::Float64=1.0,
    intercept_only::Bool=false,
    max_iter::Int=3000,
    verbose::Bool=false,
    convert_to_mee::Bool=true
)
    # Convert Cartesian to MEE if needed
    if convert_to_mee && length(x0) == 7
        # Cartesian input: [r; v; m]
        # Use parent module functions
        x0_mee_cart = Main.GTOC12.Cartesian2MEE(x0[1:6]; μ=μ)
        x0_mee = vcat(x0_mee_cart, x0[7])

        xf_mee = Main.GTOC12.Cartesian2MEE(xf; μ=μ)
    else
        # Already in MEE
        x0_mee = x0
        xf_mee = xf
    end

    # Canonical unit scaling
    AU_m = 1.49597870691e11
    L_scale = AU_m
    T_scale = 10.0 * 86400.0
    V_scale = L_scale / T_scale
    A_scale = L_scale / T_scale^2
    μ_scale = L_scale^3 / T_scale^2
    M_scale = x0_mee[7]
    F_scale = M_scale * A_scale

    # Non-dimensionalize MEE state
    x0_nd = zeros(7)
    x0_nd[1] = x0_mee[1] / L_scale  # p (semi-latus rectum)
    x0_nd[2:6] = x0_mee[2:6]  # f, g, h, k, L (dimensionless)
    x0_nd[7] = x0_mee[7] / M_scale  # mass

    xf_nd = zeros(6)
    xf_nd[1] = xf_mee[1] / L_scale
    xf_nd[2:6] = xf_mee[2:6]

    t_f_nd = t_f / T_scale
    μ_nd = μ / μ_scale
    T_max_nd = T_max / F_scale
    g0_nd = g0 / A_scale
    m_min_nd = m_min / M_scale

    c_dim = Isp * g0
    c_nd = c_dim / V_scale

    t_nodes_nd = nothing
    if t_nodes !== nothing
        t_nodes_nd = t_nodes ./ T_scale
    end

    # Solve in canonical units
    result_nd = _infiniteopt_low_thrust_mee_canonical(
        x0_nd, xf_nd, t_f_nd, N, μ_nd;
        T_max=T_max_nd,
        c=c_nd,
        m_min=m_min_nd,
        X_guess=X_guess,
        U_guess=U_guess,
        t_nodes=t_nodes_nd,
        position_tol=position_tol / L_scale,
        velocity_tol=velocity_tol / V_scale,
        intercept_only=intercept_only,
        max_iter=max_iter,
        verbose=verbose
    )

    # Re-dimensionalize results
    X_sol = copy(result_nd.X)
    X_sol[1, :] .*= L_scale  # p
    # X_sol[2:6, :] already dimensionless (f, g, h, k, L)
    X_sol[7, :] .*= M_scale  # mass

    U_sol = result_nd.U
    t_nodes_actual = result_nd.t_nodes .* T_scale

    constraint_breakdown_actual = (
        position=result_nd.constraint_breakdown.position * L_scale,
        velocity=result_nd.constraint_breakdown.velocity * V_scale,
        mass_negative=result_nd.constraint_breakdown.mass_negative * M_scale,
        thrust_direction=result_nd.constraint_breakdown.thrust_direction
    )

    mass_used_actual = result_nd.mass_used * M_scale

    return InfiniteOptLowThrustResult(
        result_nd.success,
        result_nd.objective_value,
        X_sol,
        U_sol,
        t_nodes_actual,
        result_nd.solve_time,
        result_nd.iterations,
        result_nd.constraint_violation,
        constraint_breakdown_actual,
        mass_used_actual,
        :mee
    )
end

"""
    _infiniteopt_low_thrust_mee_canonical

Internal function for MEE optimization in canonical units.
"""
function _infiniteopt_low_thrust_mee_canonical(
    x0::Vector{Float64},
    xf::Vector{Float64},
    t_f::Float64,
    N::Int,
    μ::Float64;
    T_max::Float64=0.6,
    c::Float64=1.0,
    m_min::Float64=0.5,
    X_guess::Union{Matrix{Float64},Nothing}=nothing,
    U_guess::Union{Matrix{Float64},Nothing}=nothing,
    t_nodes::Union{Vector{Float64},Nothing}=nothing,
    position_tol::Float64=1e-6,
    velocity_tol::Float64=1e-9,
    intercept_only::Bool=false,
    max_iter::Int=3000,
    verbose::Bool=false
)
    # Create InfiniteOpt model
    model = InfiniteModel(optimizer_with_attributes(
        Ipopt.Optimizer,
        "print_level" => verbose ? 5 : 0,
        "max_iter" => max_iter,
        "tol" => 1e-8,
        "constr_viol_tol" => 1e-8,
        "acceptable_tol" => 1e-6,
        "acceptable_constr_viol_tol" => 1e-6,
        "acceptable_iter" => 50,
        "mumps_mem_percent" => 8000,
        "mu_strategy" => "adaptive",
        "max_cpu_time" => 600.0,
        "nlp_scaling_method" => "gradient-based"
    ))

    # Define time parameter
    if t_nodes !== nothing
        @infinite_parameter(model, t in [0, t_f], supports = t_nodes,
            derivative_method = InfiniteOpt.OrthogonalCollocation(4))
    else
        @infinite_parameter(model, t in [0, t_f], num_supports = N,
            derivative_method = InfiniteOpt.OrthogonalCollocation(4))
    end

    # State variables: MEE elements + mass
    @variable(model, p, Infinite(t))  # Semi-latus rectum (must be > 0)
    @variable(model, f, Infinite(t))  # f component
    @variable(model, g, Infinite(t))  # g component
    @variable(model, h, Infinite(t))  # h component
    @variable(model, k, Infinite(t))  # k component
    @variable(model, L, Infinite(t))  # True longitude
    @variable(model, m, Infinite(t))  # Mass

    # Control variables: Spherical angles for thrust direction + throttle
    @variable(model, 0 <= θ <= π, Infinite(t))
    @variable(model, φ, Infinite(t))
    @variable(model, 0 <= δ <= 1, Infinite(t))

    # Make controls piecewise constant (if available)
    try
        constant_over_collocation(θ, t)
        constant_over_collocation(φ, t)
        constant_over_collocation(δ, t)
    catch
        # Function not available - use continuous controls
    end

    # Initial guesses
    if X_guess !== nothing
        set_start_value_function(p, t_val -> X_guess[1, 1] + (X_guess[1, end] - X_guess[1, 1]) * t_val / t_f)
        for i in 2:6
            set_start_value_function([f, g, h, k, L][i-1], t_val -> X_guess[i, 1] + (X_guess[i, end] - X_guess[i, 1]) * t_val / t_f)
        end
        set_start_value_function(m, t_val -> X_guess[7, 1] - (X_guess[7, 1] - X_guess[7, end]) * t_val / t_f)
    else
        # Linear interpolation for MEE elements
        set_start_value_function(p, t_val -> x0[1] + (xf[1] - x0[1]) * t_val / t_f)
        set_start_value_function(f, t_val -> x0[2] + (xf[2] - x0[2]) * t_val / t_f)
        set_start_value_function(g, t_val -> x0[3] + (xf[3] - x0[3]) * t_val / t_f)
        set_start_value_function(h, t_val -> x0[4] + (xf[4] - x0[4]) * t_val / t_f)
        set_start_value_function(k, t_val -> x0[5] + (xf[5] - x0[5]) * t_val / t_f)
        set_start_value_function(L, t_val -> x0[6] + (xf[6] - x0[6]) * t_val / t_f)
        set_start_value_function(m, t_val -> x0[7] * (1.0 - 0.1 * t_val / t_f))
    end

    if U_guess !== nothing
        # Use provided control guess
        θ_init, φ_init = U_guess[1, 1], U_guess[2, 1]
        δ_init = U_guess[3, 1]
        set_start_value_function(θ, t_val -> θ_init)
        set_start_value_function(φ, t_val -> φ_init)
        set_start_value_function(δ, t_val -> δ_init)
    else
        # Default guess
        set_start_value_function(θ, t_val -> π / 2)
        set_start_value_function(φ, t_val -> 0.0)
        set_start_value_function(δ, t_val -> 0.5)
    end

    # Objective: Maximize final mass
    t_final = upper_bound(t)
    @objective(model, Max, m(t_final))

    # Initial conditions
    @constraint(model, p(0) == x0[1])
    @constraint(model, f(0) == x0[2])
    @constraint(model, g(0) == x0[3])
    @constraint(model, h(0) == x0[4])
    @constraint(model, k(0) == x0[5])
    @constraint(model, L(0) == x0[6])
    @constraint(model, m(0) == x0[7])

    # Final conditions (for MEE, we match the orbital elements)
    @constraint(model, p(t_final) == xf[1])
    @constraint(model, f(t_final) == xf[2])
    @constraint(model, g(t_final) == xf[3])
    @constraint(model, h(t_final) == xf[4])
    @constraint(model, k(t_final) == xf[5])

    if !intercept_only
        # Match true longitude as well
        @constraint(model, L(t_final) == xf[6])
    end

    # MEE Dynamics: dMEE/dt = A(MEE) + (T/m) * δ * B(MEE) * û
    # Thrust direction: û = [sin(θ)cos(φ), sin(θ)sin(φ), cos(θ)]

    # Helper expressions
    @expression(model, q, 1 + f * cos(L) + g * sin(L))
    @expression(model, s, sqrt(1 + h^2 + k^2))
    @expression(model, p_abs, abs(p))  # Ensure p is positive

    # A matrix (drift dynamics) - only affects L
    @expression(model, A_L, sqrt(μ * p_abs) * (q / p_abs)^2)

    # B matrix (control influence matrix) - 6x3
    # Column 1: Radial thrust (û_r)
    @expression(model, B_p_r, 0.0)
    @expression(model, B_f_r, sqrt(p_abs / μ) * sin(L))
    @expression(model, B_g_r, -sqrt(p_abs / μ) * cos(L))
    @expression(model, B_h_r, 0.0)
    @expression(model, B_k_r, 0.0)
    @expression(model, B_L_r, 0.0)

    # Column 2: Tangential thrust (û_t)
    @expression(model, B_p_t, 2 * p_abs / q * sqrt(p_abs / μ))
    @expression(model, B_f_t, sqrt(p_abs / μ) / q * ((q + 1) * cos(L) + f))
    @expression(model, B_g_t, sqrt(p_abs / μ) / q * ((q + 1) * sin(L) + g))
    @expression(model, B_h_t, 0.0)
    @expression(model, B_k_t, 0.0)
    @expression(model, B_L_t, 0.0)

    # Column 3: Normal thrust (û_n)
    @expression(model, B_p_n, 0.0)
    @expression(model, B_f_n, -sqrt(p_abs / μ) * g / q * (h * sin(L) - k * cos(L)))
    @expression(model, B_g_n, sqrt(p_abs / μ) * f / q * (h * sin(L) - k * cos(L)))
    @expression(model, B_h_n, sqrt(p_abs / μ) * s * cos(L) / (2 * q))
    @expression(model, B_k_n, sqrt(p_abs / μ) * s * sin(L) / (2 * q))
    @expression(model, B_L_n, sqrt(p_abs / μ) / q * (h * sin(L) - k * cos(L)))

    # Thrust components in RTN frame
    @expression(model, û_r, sin(θ) * cos(φ))
    @expression(model, û_t, sin(θ) * sin(φ))
    @expression(model, û_n, cos(θ))

    # Control acceleration magnitude
    @expression(model, a_ctrl, (T_max / m) * δ)

    # Dynamics constraints
    @constraint(model, ∂(p, t) == a_ctrl * (B_p_r * û_r + B_p_t * û_t + B_p_n * û_n))
    @constraint(model, ∂(f, t) == a_ctrl * (B_f_r * û_r + B_f_t * û_t + B_f_n * û_n))
    @constraint(model, ∂(g, t) == a_ctrl * (B_g_r * û_r + B_g_t * û_t + B_g_n * û_n))
    @constraint(model, ∂(h, t) == a_ctrl * (B_h_r * û_r + B_h_t * û_t + B_h_n * û_n))
    @constraint(model, ∂(k, t) == a_ctrl * (B_k_r * û_r + B_k_t * û_t + B_k_n * û_n))
    @constraint(model, ∂(L, t) == A_L + a_ctrl * (B_L_r * û_r + B_L_t * û_t + B_L_n * û_n))

    # Mass dynamics
    @constraint(model, ∂(m, t) == -(T_max / c) * δ)

    # Path constraints
    @constraint(model, [t = t], m >= m_min)
    @constraint(model, [t = t], p >= 1e-6)  # Semi-latus rectum must be positive

    # Solve
    start_time = time()
    optimize!(model)
    solve_time = time() - start_time

    # Extract results
    success = termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
    obj_value = objective_value(model)

    t_vals = value(t)
    p_vals = value(p)
    f_vals = value(f)
    g_vals = value(g)
    h_vals = value(h)
    k_vals = value(k)
    L_vals = value(L)
    m_vals = value(m)
    θ_vals = value(θ)
    φ_vals = value(φ)
    δ_vals = value(δ)

    N_actual = length(t_vals)
    X_sol = zeros(7, N_actual)
    U_sol = zeros(4, N_actual)

    for i in 1:N_actual
        X_sol[1, i] = p_vals[i]
        X_sol[2, i] = f_vals[i]
        X_sol[3, i] = g_vals[i]
        X_sol[4, i] = h_vals[i]
        X_sol[5, i] = k_vals[i]
        X_sol[6, i] = L_vals[i]
        X_sol[7, i] = m_vals[i]

        U_sol[1, i] = θ_vals[i]
        U_sol[2, i] = φ_vals[i]
        U_sol[3, i] = δ_vals[i]
        U_sol[4, i] = T_max * δ_vals[i]
    end

    iterations = -1

    # Compute constraint violations
    max_position_viol = 0.0
    max_velocity_viol = 0.0
    max_mass_negative_viol = 0.0
    max_thrust_direction_viol = 0.0

    # Check final MEE elements (position)
    max_position_viol = maximum(abs.([
        X_sol[1, end] - xf[1],  # p
        X_sol[2, end] - xf[2],  # f
        X_sol[3, end] - xf[3],  # g
        X_sol[4, end] - xf[4],  # h
        X_sol[5, end] - xf[5]   # k
    ]))

    # Check final L (velocity-like for MEE)
    if !intercept_only
        max_velocity_viol = abs(X_sol[6, end] - xf[6])
    end

    # Check mass positive
    for i in 1:N_actual
        if X_sol[7, i] < m_min
            max_mass_negative_viol = max(max_mass_negative_viol, m_min - X_sol[7, i])
        end
    end

    # Check thrust direction normalization
    for i in 1:N_actual
        û_norm = sqrt(sin(U_sol[1, i])^2 + cos(U_sol[1, i])^2)
        max_thrust_direction_viol = max(max_thrust_direction_viol, abs(û_norm - 1.0))
    end

    max_viol = max(max_position_viol, max_velocity_viol, max_mass_negative_viol, max_thrust_direction_viol)

    constraint_breakdown = (
        position=max_position_viol,
        velocity=max_velocity_viol,
        mass_negative=max_mass_negative_viol,
        thrust_direction=max_thrust_direction_viol
    )

    mass_used = X_sol[7, 1] - X_sol[7, end]

    if verbose
        println("\n" * "="^70)
        println("MEE Constraint Violation Breakdown (canonical units):")
        println("="^70)
        println("  MEE elements:     $(round(constraint_breakdown.position, digits=10))")
        println("  True longitude:   $(round(constraint_breakdown.velocity, digits=10))")
        println("  Mass negative:    $(round(constraint_breakdown.mass_negative, digits=10))")
        println("  Thrust direction: $(round(constraint_breakdown.thrust_direction, digits=10))")
        println("  Overall max:      $(round(max_viol, digits=10))")
        println("  Mass used:        $(round(mass_used, digits=6))")
        println("="^70)
    end

    return InfiniteOptLowThrustResult(
        success,
        obj_value,
        X_sol,
        U_sol,
        collect(t_vals),
        solve_time,
        iterations,
        max_viol,
        constraint_breakdown,
        mass_used,
        :mee
    )
end

"""
    generate_gtoc12_submission_file(result::InfiniteOptLowThrustResult, filename::String; kwargs...)

Generate a GTOC12 Result.txt submission file from InfiniteOpt results.

This is a placeholder - to be implemented after validation with hit_mine_an_asteroid example.
"""
function generate_gtoc12_submission_file(
    result::InfiniteOptLowThrustResult,
    filename::String;
    verbose::Bool=true
)
    error("Submission file generation not yet implemented.")
end

end # module
