using SPICE
using AstroTime
using Downloads: download
using Plots
using PlanetOrbits
import PlanetOrbits: m2au, _trueanom_from_eccanom
using OrdinaryDiffEq
using LinearAlgebra

const LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
const SPK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
const PCK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc"

## Only call downloads if you don't already have the SPICE kernels
# # Download kernels
download(LSK, "naif0012.tls")
download(PCK, "gm_de440.tpc")
download(SPK, "de440.bsp")

# Load leap seconds kernel
furnsh("naif0012.tls")

# Load a planetary ephemeris kernel
furnsh("de440.bsp")

# Load a planetary properties kernel
furnsh("gm_de440.tpc")

## Convert to Ephemeris time for use with SPICE


## Asteroid 2001 TW229 
# a (semi-major axis, AU): 2.5897261
# e (eccentricity): 0.2734625
# i (inclination, deg.): 6.40734
# ω  (argument  of  pericenter,  deg.): 264.78691
# Ω  (Right  Ascension  of  the  Ascending Node, deg.): 128.34711
# M (mean anomaly at epoch 53600 MJD, deg.): 320.47955

# Compute asteroid reference time from mean anomaly
# Because the PlanetOrbits package uses reference time exclusively
meananom_to_ref_time(mean_anom) = mean_anom/(2π)


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

function Epoch_to_SPICE_ET(epoch::Epoch)
    ET_J2000 = TDBEpoch(0days, origin=:j2000)
    ET = value(seconds(AstroTime.j2000(epoch))) - value(seconds(AstroTime.j2000(ET_J2000)))
end

## Get Earth and asteroid positions at epoch
# Convert the calendar date to ephemeris seconds past J2000
# ET_J2000 = TDBEpoch(0days, origin=:j2000)
# ETs = value.(seconds.(AstroTime.j2000.(launch_dates_ET))) .- value(seconds(AstroTime.j2000(ET_J2000)))
ET_0 = launch_dates_ET[1]

## Plot Earth and asteroid positions at epoch
# Plot planet paths
# Start with a Keplerian approximation (better compatibility with
# generation of gradients via AutoDiff)

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

function keplerian2MEE(;a, e, i, Ω, ω, ν)
    """
    Converts from classical Keplerian elements to modified equinoctial elements
    """
    p = a*(1-e^2)
    f = e*cos(Ω + ω)
    g = e*sin(Ω + ω)
    h = tan(i/2)*cos(Ω)
    k = tan(i/2)*sin(Ω)
    l = Ω + ω + ν

    [p, f, g, h, k, l]

end

function MME2keplerian(;p, f, g, h, k, l)
    # See https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf

    a = p/(1 - f^2 - g^2)
    e = sqrt(f^2 + g^2)
    i = atan( 2*sqrt(h^2 + k^2), 1-h^2-k^2)
    Ω = atan( k,h )
    ω = atan( gh-fk, fh+gk )
    ν = l - (Ω + ω)

    

    [a, e, i, Ω, ω, ν]
end

get_q(;f, g, L) = 1 + f*cos(L) + g*sin(L)
get_s(;h, k) = sqrt(1 + h^2 + k^2)
get_α(;h, k) = h^2 - k^2
get_r(;p, w) = p/w
get_w(;f::Number, g::Number, L::Number) = 1+f*cos(L)+g*sin(L)

function A_equinoctial(MEE; μ)
    # Eq 2
    p, f, g, h, k, l = MEE
    q = get_q(f=f, g=g, L=l)

    A = [0;0;0;0;0;sqrt(μ*p)*(q/p)^2]
end

function B_equinoctial(MEE; μ)
    # Eq 3
    p, f, g, h, k, l = MEE
    q = get_q(f=f, g=g, L=l)
    s = get_s(h=h, k=k)
    
    B = [               0                2*p/q*sqrt(p/μ)                                    0
         sqrt(p/μ)*sin(l)  sqrt(p/μ)*q*((q+1)*cos(l) + f) -sqrt(p/μ)*g/q*(h*sin(l) - k*cos(l))
        -sqrt(p/μ)*cos(l)  sqrt(p/μ)*q*((q+1)*sin(l) + g) -sqrt(p/μ)*g/q*(h*sin(l) - k*cos(l))
                        0                               0             sqrt(p/μ)*s*cos(l)/(2*q)
                        0                               0             sqrt(p/μ)*s*sin(l)/(2*q)
                        0                               0  sqrt(p/μ)*1/q*(h*sin(l) - k*cos(l))]
end

function MEE2Cartesian(MEE; μ)
    # Source: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
    # Equations 3a & 3b
    p, f, g, h, k, l = MEE
    
    α = get_α(h=h, k=k)
    s = get_s(h=h, k=k)
    w = get_w(;f=g, g=g, L=l) 

    r = p/w
    r⃗ = r/(s^2) * [cos(l) + α^2 * cos(l) + 2*h*k*sin(l)
                    sin(l) - α^2 * sin(l) + 2*h*k*cos(l)
                    2*(h*sin(l) - k*cos(l))]

    v⃗ = 1/s^2 * sqrt(μ/p) *[-(sin(l) +α^2*sin(l) - 2*h*k*cos(l) + g -2*f*h*k + α^2*g)
                            -(-cos(l) + α^2*sin(l) + 2*h*k*sin(l) - f +2*g*h*k + α^2*g)
                            2*(h*cos(l) + k*sin(l) + f*h + g*k)]

    vcat(r⃗, v⃗)
end

function get_control(MEE; params=p) # Get control thrust direction and magnitude [0, 1]
    μ, c, T = params # problem parameters

    x⃗ = MEE2Cartesian(MEE; μ)
    v⃗ = x⃗[4:6]

    # Just roll with tangential firing to sanity check
    û = v⃗/norm(v⃗)
    δ = 1
    
    û, δ 
end

function EOM_MEE!(ẋ, x, p, t)
    μ, c, T = p # problem parameters
    MEE = x[1:6]
    m   = x[7]

    A = A_equinoctial(MEE; μ=μ)
    B = B_equinoctial(MEE; μ=μ)
    û, δ  = get_control(MEE; params=p) # Get control thrust direction and magnitude [0, 1]

    MEE_dot = A + T*δ/m * B * û
    ṁ = -δ*T/(c)

    ẋ[:] = vcat( MEE_dot, ṁ )

end

# Start spacecraft at Earth and propagate with the tangential guidance
tspan = (0, 3600)

M_earth = PlanetOrbits.trueanom
MEE₀ = keplerian2MEE(;
a = earth_orbit.a /m2au,
e = earth_orbit.e,
i = earth_orbit.i,
Ω = earth_orbit.Ω, 
ω = earth_orbit.ω,
ν = trueanom(earth_orbit, tspan[1]) )

# Problem parameters
μ_☉ = bodvrd("Sun", "GM")[1] # Sun central body
c = Isp * g₀ # exhaust velocity


parameters = ( μ_☉, c, T_max)
prob = ODEProblem(EOM_MEE!, vcat(MEE₀, m₀), tspan, parameters)
sol = solve(prob, Tsit5())

MEE_out = [sol[1:6,i] for i = 1:length(sol)]

# Convert the Mean Equinoctial Elements to Cartesian
x_spacecraft = MEE2Cartesian.(MEE_out; μ=μ_☉)
r_spacecraft = getindex.(x_spacecraft', 1:3)'.*m2au

plot(r_spacecraft[:,1], r_spacecraft[:,2])
