export optimize_impulsive_launch

## Problem definition
# TODO: Make a type that has passthrough calls to the methods already in Altro.jl
function optimize_impulsive_launch(x₀, x_target, Δt; ΔV₀=zeros(3), V∞_max=6000.0, ET_launch=GTOC12.ET₀, asteroid_ID)
    # Problem scaling via canonical units (this makes the solver not hate you)
    x₀, CDU, CTU, μ_canonical = get_canonical_state(x₀; μ=GTOC12.μ_☉)
    canonical_state!(x_target; CDU=CDU, CTU=CTU)
    V∞_max = canonical_vel(V∞_max; CDU=CDU, CTU=CTU)
    canonical_vel!(ΔV₀; CDU=CDU, CTU=CTU)
    Δt = canonical_time(Δt; CTU)

    # Model and discretization
    target_asteroid = Asteroid(asteroid_ID)
    model = ImpulsiveSpacecraft(ET_launch=ET_launch, target_asteroid=target_asteroid, μ=μ_canonical)
    n, m = RD.dims(model)
    tf = Δt  # final time (sec)
    N = 2          # number of knot points, 2 for simple start/end impulses
    dt = tf / (N - 1)  # time step (sec)

    # Objective
    # Note: May have to tune on a segment-by-segment basis
    x0 = SVector{length(x₀)}(x₀)# initial state
    xf = SVector{length(x_target)}(x_target) # final state

    # Set up costs.
    # TODO: Convert costs to just Σ(ΔV)
    Q = Diagonal(@SVector zeros(n))         # Don't penalize intermediate states at all.
    R = 1e-10 * Diagonal(@SVector ones(m))  # Don't Penalize ΔV (control) since we get it for free from launch vehicle
    Qf = Diagonal(@SVector [0; 0; 0; 1; 1; 1])  # Only penalize final velocity (treat ΔV to match at the end as part of objective)
    obj = LQRObjective(Q, R, Qf, xf, N)

    # Add constraints
    conSet = ConstraintList(n, m, N)
    add_constraint!(conSet, GoalConstraint(xf, 1:3), N) # Just constrain the position at the end
    bnd = BoundConstraint(n, m, u_min=-V∞_max / sqrt(m) * ones(m), u_max=V∞_max / sqrt(m) * ones(m)) # 6 km/s (RSS) limit on hyperbolic excess velocity
    add_constraint!(conSet, bnd, 1:N-1)

    # Initial guess
    u0 = SVector{length(ΔV₀)}(ΔV₀) # Try warm-starting with shooting method solution?
    # u0 = @SVector zeros(m)

    # Set up problem
    prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)
    initial_controls!(prob, u0)
    rollout!(RD.InPlace(), prob)


    # Uh, solve the problem?
    # See options here: https://github.com/RoboticExplorationLab/Altro.jl#list-of-options
    opts = SolverOptions(;
        penalty_initial=1e6,
        penalty_scaling=1e6,
        constraint_tolerance=1e6,
        cost_tolerance=1e6
    )

    solver = ALTROSolver(prob, opts)
    solve!(solver)

    # Extract states and controls
    states_out = states(prob)
    controls_out = controls(prob)
    times_out = gettimes(prob)

    # Redimensionalize input states InPlace
    redimensionalize_state!(x_target; CDU=CDU, CTU=CTU)
    redimensionalize_vel!(ΔV₀; CDU=CDU, CTU=CTU)

    # Redimensionalize output based on canonical units
    redimensionalize_state!.(states_out; CDU=CDU, CTU=CTU)
    redimensionalize_vel!.(controls_out; CDU=CDU, CTU=CTU)
    times_out = redimensionalize_time.(times_out; CTU=CTU)

    return states_out, controls_out, times_out
end