using GTOC1
using DifferentialEquations, Plots

## Furnish relevant SPICE kernels
furnish_all_kernels()

## Asteroid 2001 TW229 
# a (semi-major axis, AU): 2.5897261
# e (eccentricity): 0.2734625
# i (inclination, deg.): 6.40734
# ω  (argument  of  pericenter,  deg.): 264.78691
# Ω  (Right  Ascension  of  the  Ascending Node, deg.): 128.34711
# M (mean anomaly at epoch 53600 MJD, deg.): 320.47955

asteroid_orbit = orbit(
    a=2.5897261,          # semi major axis (AU)
    M=1.0,                # primary mass (solar masses)
    e=0.2734625,          # eccentricity
    ω=deg2rad(264.78691), # argument of periapsis (radians)
    i=deg2rad(6.40734),   # inclination (radians)
    Ω=deg2rad(128.34711), # inclination (radians)
    τ=meananom_to_ref_time( deg2rad(320.47955) ),                 # Reference time since perihelion
    tref=53600            # Modified Julian Date
)


## Spacecraft
m₀ = 1500; # kg
Isp = 2500; # seconds
g₀ = 9.80665 # Wikipedia
T_max = .04; # Newtons

## Mission
launch_dates_ET = TDBEpoch.([3653days; 10958days], origin=:j2000); # TDB is effectively ET
ToF_max = 30*365*24*3600; # seconds
min_radius_AU = 0.2; # Minimum solar radius allowed
v∞_launch = 2.5 * 1000; # m/s, launch velocity relative to Earth
# AU=1.4959787066e+008  km
# g0=9.80665 m/s2

## Flyby constraints
#                                            Mercury        Venus          Earth             Mars          Jupiter           Saturn         Sun
#  Gravitational Constant, km^3/sec^2         22321         324860     398601.19            42828.3       126700000        37900000     1.32712428e+011 
#  Min peri radius allowed during fly-by, km    2740           6351           6678              3689           600000           70000           N/A   

# Objective function:
# J = m_f*| U_rel ⋅ v_asteroid |


## Get Earth and asteroid positions at epoch
# Convert the calendar date to ephemeris seconds past J2000
# ET_J2000 = TDBEpoch(0days, origin=:j2000)
# ETs = value.(seconds.(AstroTime.j2000.(launch_dates_ET))) .- value(seconds(AstroTime.j2000(ET_J2000)))
ET_0 = launch_dates_ET[1]

## Plot Earth and asteroid positions at epoch
# Plot planet paths
# Start with a Keplerian approximation (better compatibility with
# generation of gradients via AutoDiff)

earth_orbit = get_planet_orbit(;planet="Earth", ET=ET_0)

plot(earth_orbit,label="Earth")
plot!(asteroid_orbit, label="Asteroid")


## Just have the spacecraft fire tangentially and plot
# Couple options for dynamics:
#   1. Cartesian (with inverse square law)
#   2. Orbital elements (Map effects onto OE directly):
#       - This feels better
#       - May not be better for targetting
#       - ^May not matter for AutoDiff, that would be cool

# Let's use the Mean Equinoctial Elements:
# https://ai.jpl.nasa.gov/public/documents/papers/AAS-22-015-Paper.pdf
# x= [p,f,g,h,k,l]


# Start spacecraft at Earth and propagate with the tangential guidance
tspan = (0, 365*24*3600)

M_earth = PlanetOrbits.trueanom
MEE₀ = keplerian2MEE(;
a = earth_orbit.a /m2au,
e = earth_orbit.e,
i = earth_orbit.i,
Ω = earth_orbit.Ω, 
ω = earth_orbit.ω,
ν = trueanom(earth_orbit, tspan[1]) ) # Start off at Earth

# Problem parameters
μ_☉ = bodvrd("Sun", "GM")[1] * (1000)^3 # Sun central body, km³/s² → m³/s²
c = Isp * g₀ # exhaust velocity

struct GTOCProblemParams
    μ::Float64      # Central body (Sun) standard gravitational parameter [m³/s²]
    c::Float64      # Thruster exhaust velocity [m/s]
    T_max::Float64  # Max thruster force [N]
    ΔV_LV_inrt::Vector{Float64} # Initial ΔV w.r.t. Earth from Launch Vehicle (not to exceed problem value) [m/s]

    # kwarg constructor
    function GTOCProblemParams(;μ, c, T_max, ΔV_LV_inrt)
        new(μ, c, T_max, ΔV_LV_inrt)
    end
end

# Initial DV kick
DV₀ = [0; v∞_launch; 0]

parameters = GTOCProblemParams( μ=μ_☉, c=c, T_max=T_max, ΔV_LV_inrt=DV₀)
prob = ODEProblem(EOM_MEE!, vcat(MEE₀, m₀), tspan, parameters, callback=LV_callback)
sol = solve(prob, tstops =1, alg_hints = [:stiff], reltol = 1e-10, abstol = 1e-6)

MEE_out = [sol[1:6,i] for i = 1:lastindex(sol)]

# Convert the Mean Equinoctial Elements to Cartesian
x_spacecraft = MEE2Cartesian.(MEE_out; μ=μ_☉)
r_spacecraft = getindex.(x_spacecraft', 1:3)'.*m2au

plot!(r_spacecraft[:,1], r_spacecraft[:,2], label="Spacecraft")
xflip!(false)
