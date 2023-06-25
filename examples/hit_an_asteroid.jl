using GTOC12
using Plots
using SPICE
using ForwardDiff
using LinearAlgebra

# Import asteroids
asteroid_df = get_asteroid_df()

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(asteroid_df.sma .- 1))

function get_asteroid_state(asteroid, ET)
    # Change phasing to be time-based
    n = mean_motion(;a=asteroid.sma, μ=μ_☉) # mean motion

    # Propagate mean anomaly to ET
    M = asteroid.mean_anom + n * (ET - GTOC12.ET₀)

    # Calculate true anomaly
    E = M2EH(; M=M, ecc=asteroid.ecc)
    ν = EH2ν(; E_or_H=E, ecc=asteroid.ecc)

    # Get asteroid position at all true anomalies
    x_asteroid = COE2RV(; COE=
        [asteroid.sma,
            asteroid.ecc,
            asteroid.inc,
            asteroid.LAN,
            asteroid.argperi,
            ν],
        μ_CB_or_CB_name=μ_☉)
end

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

# Use that asteroids orbital elements and sweep through fast variable to plot
# TODO: Move to own plotting file
function plot_asteroid_from_df_row(asteroid; ET_in=[GTOC12.ET₀, GTOC12.ET₀ + 365 * 24 * 3600])

    # Accept start/end time or vector of times
    ETs = vec_from_span(ET_in)

    # Get asteroid position at all true anomalies
    x_asteroid = hcat([get_asteroid_state(asteroid, ET) for ET in ETs]...)

    plot(x_asteroid[1, :], x_asteroid[2, :], label="Asteroid" * string(asteroid.ID), aspect_ratio=:equal, color=colormap("Reds"))


end

# Plot the asteroid trajectory
asteroid = asteroid_df[ID_min, :]
plot_asteroid_from_df_row(asteroid)

# Wrap SPICE calls
function get_planet_state(planet, ET)
    # Use SPICE to get position
    ref = "ECLIPJ2000"
    abcorr = "NONE"
    obs = GTOC12.default_CB_str

    # TODO: Use SPICE Int IDs instead of string
    r_planet = spkezr(uppercase(planet), ET, ref, abcorr, uppercase(obs))[1]
    r_planet *= 1000 # km → m
end

# Plot Earth trajectory
planet = "EARTH"
function plot_planet!(; planet="EARTH", ET_in=[GTOC12.ET₀, GTOC12.ET₀ + 365 * 24 * 3600])

    # Accept start/end time or vector of times
    ETs = vec_from_span(ET_in)

    x_planet = hcat([get_planet_state(planet, ET) for ET in ETs]...)

    plot!(x_planet[1, :], x_planet[2, :], label=planet, color=colormap("Blues"))

end

# Plot Earth
plot_planet!(planet="EARTH")


# Try out shooting method to hit asteroid
DV₀ = [-GTOC12.v∞_max; 0; 0]
x₀ = get_planet_state("EARTH", GTOC12.ET₀) + [0;0;0;DV₀[:]]

t_shoot = 1 / 2 * 365 * 24 * 3600

# Fixed-time shooting
# 1. Get asteroid position (target)
ET_target = GTOC12.ET₀ + t_shoot
x_target = get_asteroid_state(asteroid, ET_target)

# 2. Propagate spacecraft from initial state
xₜ = propagate_keplerian(x₀, t_shoot)

# 3. Get STM: use ForwardDiff to get dr_f/dv_0
# Actual "cost" is miss distance
# Need f(v₀) = dr
function calculate_asteroid_miss_distance(x_aug)
    # Initial guess state values
    x₀ = x_aug[1:6]
    t = x_aug[7]
    rₜ_des = x_aug[8:10]

    # TODO: Augment target state to calculate sensitivity w.r.t. initial ET

    # ET_target = GTOC12.ET₀ + t
    # asteroid = get_asteroid_df()[ID, :]
    # x_target = get_asteroid_state(asteroid, ET_target)
    xₜ = propagate_keplerian(x₀, t)

    cost = xₜ[1:3] - rₜ_des
end

# TODO: Think about clever way to package problem into a type
# such that can have one giant vector for calculating adjoint_sensitivities
# 
# 5. Iterate until convergence
MAX_ITER = 50
num_iter = 0
drₜ = x₀[1:3]
while num_iter < MAX_ITER && norm(drₜ) > GTOC12.POS_ABS_TOL
    # Form vector used by ForwardDiff problem
    x_aug = vcat(x₀, t_shoot, x_target[1:3])

    # Calculate miss distance from current propagation
    global drₜ = calculate_asteroid_miss_distance(x_aug)
    println(drₜ)

    # Calculate sensitivity of solution to all parameters
    dr_dx_aug = ForwardDiff.jacobian(calculate_asteroid_miss_distance, x_aug)

    # Get partials w.r.t. design variables
    # (In this case, initial velocity and time of flight)
    ∂rₜ_∂v₀ = dr_dx_aug[:, 4:6] # Variation of final position w.r.t. initial velocity
    ∂rₜ_∂TOF = dr_dx_aug[:, 7]  # Variation of final position w.r.t. time-of-flight

    # ^^^TODO: handle indices more intelligently for packing/unpacking


    # Update initial state via Newton's method
    # See Pavlak Thesis Eq. 3.12
    dv₀ = -∂rₜ_∂v₀\drₜ
    dTOF = -∂rₜ_∂TOF\drₜ

    global x₀[4:6] += dv₀

    # Variable time stuff isn't working interestingly
    # global t_shoot += dTOF

    global num_iter += 1

end

# Plot spacecraft trajectory
function plot_coast!(x₀, tspan)
    tvec = vec_from_span(tspan)
    states = hcat([propagate_keplerian(x₀, t) for t in tvec]...)
    
    plot!(states[1, :], states[2, :], label="coast", color=colormap("Greens"))

end

plot_coast!(x₀, t_shoot)
