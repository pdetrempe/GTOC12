export Planet, Venus, Earth, Mars

struct Planet <: CelestialBody
    # Planet name
    name::String

    # Planet Orbital elements at epoch
    sma::Float64
    ecc::Float64
    inc::Float64
    LAN::Float64
    argperi::Float64
    mean_anom::Float64
    oe::SVector{6,Float64}

    # Planet epoch
    ET₀::Float64

    # Cartesian state
    x_ET₀::SVector{6,Float64}

    # Planet/flyby parameters
    μ::Float64
    r_peri_min::Float64

    # Calculate relevant quantities at construction
    function Planet(name::String, sma::Float64, ecc::Float64, inc::Float64, LAN::Float64, argperi::Float64, mean_anom::Float64, ET₀::Float64, μ::Float64, r_peri_min::Float64)
        ν = M2ν(; M=mean_anom, ecc=ecc)
        oe = @SVector [sma, ecc, inc, LAN, argperi, ν]
        x_ET₀ = COE2RV( COE=oe )
        return new( name, sma, ecc, inc, LAN, argperi, mean_anom, oe, ET₀, x_ET₀, μ, r_peri_min)
    end
end

function Planet(;name, sma, ecc, inc, LAN, argperi, mean_anom, ET₀=GTOC12.ET₀, μ, r_peri_min )
    return Planet( name, sma, ecc, inc, LAN, argperi, mean_anom, ET₀, μ, r_peri_min )
end

# Initialize all planets
const Venus = Planet(name="Venus", 
sma = 1.08208010521e8 * 1000,  # Convert to m
ecc = 6.72988099539e-3, 
inc = deg2rad( 3.39439096544 ), 
LAN = deg2rad( 7.65796397775e1 ),
argperi = deg2rad( 5.51107191497e1 ), 
mean_anom = deg2rad( 1.11218416921e1 ), 
ET₀ = GTOC12.ET₀,
μ = 3.24858592000e5 * 1000^3, # Convert to m^3/s^2
r_peri_min = 6351.0 *1000# Convert to m
)

const Earth = Planet(name="Earth", 
sma = 1.49579151285e8 * 1000,  # Convert to m
ecc = 1.65519129162e-2, 
inc = deg2rad( 4.64389155500e-3  ), 
LAN = deg2rad( 1.98956406477e2 ),
argperi = deg2rad( 2.62960364700e2 ), 
mean_anom = deg2rad( 3.58039899470e2  ), 
ET₀ = GTOC12.ET₀,
μ = 3.98600435436e5 * 1000^3, # Convert to m^3/s^2
r_peri_min = 6678.0*1000 # Convert to m
)

const Mars = Planet(name="Mars", 
sma = 12.27951663551e8 * 1000,  # Convert to m
ecc = 9.33662184095e-2, 
inc = deg2rad( 1.84693231241  ), 
LAN = deg2rad( 4.94553142513e1 ),
argperi = deg2rad( 2.86731029267e2 ), 
mean_anom = deg2rad( 2.38232037154e2  ), 
ET₀ = GTOC12.ET₀,
μ = 4.28283752140e4 * 1000^3, # Convert to m^3/s^2
r_peri_min = 3689.0 *1000 # Convert to m
)

