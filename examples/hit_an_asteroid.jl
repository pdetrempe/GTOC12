using GTOC12
using Plots

# Import asteroids
asteroid_df = get_asteroid_df()

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(asteroid_df.sma .- 1))

# Use that asteroids orbital elements and sweep through fast variable to plot
# TODO: Move to own plotting file
function plot_asteroid_from_df(asteroid)

    TA_vec = collect(range(0, 2π, 100))

    # Get asteroid position at all positions
    # TODO: vectorize conversions

    x_asteroid = hcat([COE2RV(; COE=
        [asteroid.sma,
            asteroid.ecc,
            asteroid.inc,
            asteroid.LAN,
            asteroid.argperi,
            ν],
        μ_CB_or_CB_name=μ_☉) for ν in TA_vec]...)

    plot(x_asteroid[1, :], x_asteroid[2, :], label="Asteroid" * string(asteroid.ID), aspect_ratio=:equal)


end

# Plot the asteroid trajectory
asteroid = asteroid_df[ID_min, :]
plot_asteroid_from_df(asteroid)

# Plot Earth trajectory


