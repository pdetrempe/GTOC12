using Plots

export plot_coast, plot_coast!, plot_body, plot_body!

# Plot spacecraft trajectory
function plot_coast!(x₀, tspan; kwargs...)
    tvec = vec_from_span(tspan)
    states = hcat([propagate_universal(x₀, t) for t in tvec]...)

    plot!(states[1, :], states[2, :]; kwargs...)

end

function plot_coast(x₀, tspan; kwargs...)
    tvec = vec_from_span(tspan)
    states = hcat([propagate_universal(x₀, t) for t in tvec]...)

    plot(states[1, :], states[2, :]; kwargs...)

end

# Use that asteroids orbital elements and sweep through fast variable to plot
function get_body_label(body::CelestialBody)
    if typeof(body) == Asteroid
        ID = body.ID
        label = "Asteroid $ID"
    elseif typeof(body) == Planet # is a planet
        label = body.name
    else
        label = "Body"
    end
end

function plot_body(body::CelestialBody; ETs=[GTOC12.ET₀, GTOC12.ET₀ + 365 * 24 * 3600], kwargs...)
    x_body = get_states_from_ETs(body; ETs=vec_from_span(ETs))
    label = get_body_label(body)
    plot(x_body[1, :], x_body[2, :]; label=label, aspect_ratio=:equal, color=colormap("Reds"), kwargs...)
end

function plot_body!(body::CelestialBody; ETs=[GTOC12.ET₀, GTOC12.ET₀ + 365 * 24 * 3600], kwargs...)
    x_body = get_states_from_ETs(body; ETs=vec_from_span(ETs))
    label = get_body_label(body)
    plot!(x_body[1, :], x_body[2, :]; label=label, aspect_ratio=:equal, color=colormap("Reds"), kwargs...)
end

