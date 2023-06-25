using Plots

export plot_coast!, plot_asteroid_from_df_row, plot_planet!

# Plot spacecraft trajectory
function plot_coast!(x₀, tspan; kwargs...)
    tvec = vec_from_span(tspan)
    states = hcat([propagate_keplerian(x₀, t) for t in tvec]...)
    
    plot!(states[1, :], states[2, :]; kwargs...)

end

# Use that asteroids orbital elements and sweep through fast variable to plot
function plot_asteroid_from_df_row(asteroid; ET_in=[GTOC12.ET₀, GTOC12.ET₀ + 365 * 24 * 3600], kwargs...)

    # Accept start/end time or vector of times
    ETs = vec_from_span(ET_in)

    # Get asteroid position at all true anomalies
    x_asteroid = hcat([get_asteroid_state(asteroid, ET) for ET in ETs]...)

    plot(x_asteroid[1, :], x_asteroid[2, :]; label="Asteroid" * string(asteroid.ID), aspect_ratio=:equal, color=colormap("Reds"), kwargs...)


end

# Plot planet
function plot_planet!(; planet="EARTH", ET_in=[GTOC12.ET₀, GTOC12.ET₀ + 365 * 24 * 3600], kwargs...)

    # Accept start/end time or vector of times
    ETs = vec_from_span(ET_in)

    x_planet = hcat([get_planet_state(planet, ET) for ET in ETs]...)

    plot!(x_planet[1, :], x_planet[2, :]; kwargs...)

end
