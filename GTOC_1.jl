using SPICE
using AstroTime
using Downloads: download
using Plots
using PlanetOrbits
import PlanetOrbits: m2au

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
ETs = Epoch_to_SPICE_ET.(launch_dates_ET)
ET_0 = launch_dates_ET[1]

# Get the position of Mars at `et` w.r.t. Earth
spkpos_out =  spkpos.("Earth", ETs, "ECLIPJ2000", "none", "Sun")

Earth_pos = [x[1] for x in spkpos_out]
Earth_pos_array = permutedims(hcat(Earth_pos...))

## Plot Earth and asteroid positions at epoch
# Plot planet paths
# Start with a Keplerian approximation (better compatibility with
# generation of gradients via AutoDiff)

function get_planet_orbit(;planet::String, ET::Epoch, frame="ECLIPJ2000", CB="Sun")
(state, _) = spkezr(planet, Epoch_to_SPICE_ET(ET), frame, "none", CB)
μ_CB = bodvrd(CB, "GM")[1]
μ_☉ = bodvrd("Sun", "GM")[1]
OE = oscelt(state, Epoch_to_SPICE_ET(ET), μ_CB)

e = OE[2]
a = OE[1]/(1-e) * 1000 # convert from km to m
i = OE[3]
Ω = OE[4]
ω = OE[5]
M = OE[6]


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

plot(Earth_pos_array[:,1], Earth_pos_array[:,2], Marker='o')

