export meananom_to_ref_time, get_osc_elt, get_planet_orbit

# Compute asteroid reference time from mean anomaly
# Because the PlanetOrbits package uses reference time exclusively
meananom_to_ref_time(mean_anom) = mean_anom/(2π)

# Get osculating elements for a planet
function get_osc_elt(;planet::String, ET::Epoch, frame="ECLIPJ2000", CB="Sun")
    (state, _) = spkezr(planet, Epoch_to_SPICE_ET(ET), frame, "none", CB)
    μ_CB = bodvrd(CB, "GM")[1]
    OE = oscelt(state, Epoch_to_SPICE_ET(ET), μ_CB)
end


function get_planet_orbit(;planet::String, ET::Epoch, frame="ECLIPJ2000", CB="Sun")
    OE = get_osc_elt(;planet=planet, ET=ET, frame=frame, CB=CB)

    e = OE[2]
    a = OE[1]/(1-e) * 1000 # convert from km to m
    i = OE[3]
    Ω = OE[4]
    ω = OE[5]
    M = OE[6]

    μ_CB = bodvrd(CB, "GM")[1]
    μ_☉ = bodvrd("Sun", "GM")[1]


    orbit(    
    a=a * m2au,            # semi major axis (AU)
    M=μ_CB/μ_☉,                   # primary mass (solar masses)
    e=e,                          # eccentricity
    ω=ω,                          # argument of periapsis (radians)
    i=i,                          # inclination (radians)
    Ω=Ω,                          # RAAN (radians)
    τ=meananom_to_ref_time( M ),  # Reference time since perihelion
    tref=value(modified_julian(ET)))

end

