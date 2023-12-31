using SPICE

export μ_☉, au2m, ET₀

# From Problem statement appendix
const μ_☉ = 1.32712440018e11 * (1000)^3 # Sun central body, km³/s² → m³/s²
const au2m = 1.49597870691e11
const v∞_max = 6 * 1000; # m/s, hyperbolic excess velocity relative to Earth
const POS_ABS_TOL = 1e6; # 1000 km miss distance tolerance for asteroid events
const VEL_ABS_TOL = 1;   # 1 m/s tolerance
const MASS_ABS_TOL = .001; # .001 kg tolerance
const g0 = 9.80665 # m/s^2
const Isp = 4000   # seconds 
const T_max = 0.6 # Max Thrust 
const TT_start = "2035, Jan 1, 00:00:00.0000 (TDT)" # TT is TDT in SPICE
const MJD_0 = 64328.0
const year2day = 365.25
const day2sec =  24.0 * 3600
const year2sec = year2day * day2sec