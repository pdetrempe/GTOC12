using GTOC12
using Plots

# Import asteroids
asteroid_df = get_asteroid_df()

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(asteroid_df.sma .- 1))

# Use that asteroids orbital elements and sweep through fast variable to plot
function plot_asteroid_from_df(asteroid)

    TA_vec = collect(range(0, 2π, 100))

    # Get asteroid position at all positions
    # TODO: vectorize conversions

    x_asteroid = hcat([COE2RV(; COE=[asteroid.sma * au2m,
            asteroid.ecc,
            deg2rad(asteroid.inc),
            deg2rad(asteroid.LAN),
            deg2rad(asteroid.argperi),
            ν], μ_CB_or_CB_name=μ_☉) for ν in TA_vec]...)

    plot(x_asteroid[1,:], x_asteroid[2,:], label="Asteroid" * string(asteroid.ID))


end

# Plot the asteroid trajectory
plot_asteroid_from_df(asteroid)
