export vec_from_span, get_canonical_position, get_canonical_velocity,  get_canonical_state, 
canonical_time, canonical_vel, canonical_pos!, canonical_vel!, canonical_state!,
redimensionalize_state!, redimensionalize_pos!, redimensionalize_vel!, redimensionalize_time

# Helper function to add filler points to start/end span
function vec_from_span(span; num_points=100)
    if length(span) == 1
        vec_out = collect(range( 0, span[1], num_points))

    elseif length(span) == 2
        vec_out = collect(range( span[1], span[2], num_points))
        
    else
        return span
    end
end

# Helper function to get states for plotting
function get_states_from_ETs(body::CelestialBody; ETs)
    x_planet = hcat([get_body_state(body; ET=ET) for ET in ETs]...)
end

#---------------------- Utilities for converting to/from canonical states
get_canonical_position(r; CDU) = r/CDU
get_canonical_velocity(v; CDU, CTU) = v/(CDU/CTU)
get_canonical_time(t; CTU) = t/CTU


function get_canonical_state(x_dimensional, CDU, CTU)
    r_dim = x_dimensional[1:3]
    v_dim = x_dimensional[4:6]
    [get_canonical_position(r_dim; CDU=CDU) ; get_canonical_velocity(v_dim; CDU=CDU, CTU=CTU)] # Non-dimensionalize position/velocity values
end

function get_canonical_state(x_dimensional, CDU; μ=GTOC12.μ_☉)
    CTU = √(CDU^3/μ)   # Canonical Time Unit
    get_canonical_state(x_dimensional, CDU, CTU)
end

function get_canonical_state(x_dimensional; μ=GTOC12.μ_☉)
    r_dim = x_dimensional[1:3]
    CDU = norm(r_dim)  # Canonical Distance Unit
    CTU = √(CDU^3/μ)   # Canonical Time Unit
    μ_canonical = 1.0

    x_canonical = get_canonical_state(x_dimensional, CDU, CTU)
    return x_canonical, CDU, CTU, μ_canonical
end

# Functions for making canonical in place
function canonical_state!(x_in; CDU, CTU)
    length(x_in) == 6 || throw(ArgumentError("Invalid length of (x_in = $x_in), should be 6"))
    x_in[:] = get_canonical_state(x_in, CDU, CTU)
    nothing
end

function canonical_pos!(r_in; CDU)
    r_in .= get_canonical_position(r_in; CDU=CDU)
    nothing
end

function canonical_vel!(v_in; CDU, CTU)
    v_in .= get_canonical_velocity(v_in; CDU=CDU, CTU=CTU)
    nothing
end

function canonical_vel(v_in; CDU, CTU)
    get_canonical_velocity(v_in; CDU=CDU, CTU=CTU)
end

function canonical_time(t_in; CTU)
    get_canonical_time(t_in; CTU=CTU)
end

# Redimensionalize values by inverting canonical values
function redimensionalize_state!(x_in; CDU, CTU)
    CDU_dim = 1/CDU
    CTU_dim = 1/CTU
    canonical_state!(x_in; CDU=CDU_dim, CTU=CTU_dim)
end

function redimensionalize_pos!(r_in; CDU)
    CDU_dim = 1/CDU
    canonical_pos!(r_in; CDU=CDU_dim)
end

function redimensionalize_vel!(v_in; CDU, CTU)
    CDU_dim = 1/CDU
    CTU_dim = 1/CTU
    canonical_vel!(v_in; CDU=CDU_dim, CTU=CTU_dim)
end

function redimensionalize_time(t_in; CTU)
    CTU_dim = 1/CTU
    t_in = canonical_time(t_in; CTU=CTU_dim)
end