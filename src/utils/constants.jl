using SPICE

export μ_☉, au2m

# From Problem statement appendix
const μ_☉ = 1.32712440018e11 * (1000)^3 # Sun central body, km³/s² → m³/s²
const au2m = 1.49597870691e11
const v∞_max = 6 * 1000; # m/s, hyperbolic excess velocity relative to Earth

asteroid_df = get_asteroid_df()
const ET₀ = asteroid_df[1, :ET]