export CelestialBody, get_body_state

abstract type CelestialBody end # Type for querying celestial body states

# Function for propagating Planet state to epoch
function get_body_state(body::CelestialBody; ET)
    Δt = ET - body.ET₀
    propagate_universal(body.x_ET₀, Δt)
end