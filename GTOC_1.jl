using SPICE
using AstroTime
using Downloads: download
using Plots

const LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
const SPK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"

# Download kernels
download(LSK, "naif0012.tls")
download(SPK, "de440.bsp")

# Load leap seconds kernel
furnsh("naif0012.tls")

# Load a planetary ephemeris kernel
furnsh("de440.bsp")

## Convert to Ephemeris time for use with SPICE


## Asteroid 2001 TW229 
# a (semi-major axis, AU): 2.5897261
# e (eccentricity): 0.2734625
# i (inclination, deg.): 6.40734
# ω  (argument  of  pericenter,  deg.): 264.78691
# Ω  (Right  Ascension  of  the  Ascending Node, deg.): 128.34711
# M (mean anomaly at epoch 53600 MJD, deg.): 320.47955

## Spacecraft
m₀ = 1500; # kg
Isp = 2500; # seconds
T_max = .04; # Newtons

## Mission
launch_dates_ET = TDBEpoch.([3653days; 10958days], origin=:j2000); # TDB is effectively ET
ToF_max = 30*365*24*3600; # seconds
min_radius_AU = 0.2; # Minimum solar radius allowed
v_inf_launch = 2.5 * 1000; # m/s, launch velocity relative to Earth
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
ET_0 = TDBEpoch(0days, origin=:j2000)
ETs = value.(seconds.(AstroTime.j2000.(launch_dates_ET))) .- value(seconds(AstroTime.j2000(ET_0)))

# Get the position of Mars at `et` w.r.t. Earth
spkpos_out =  spkpos.("Earth", ETs, "ECLIPJ2000", "none", "Sun")

Earth_pos = [x[1] for x in spkpos_out]
Earth_pos_array = permutedims(hcat(Earth_pos...))

## Plot Earth and asteroid positions at epoch
plot(Earth_pos_array[:,1], Earth_pos_array[:,2], Marker='o')

