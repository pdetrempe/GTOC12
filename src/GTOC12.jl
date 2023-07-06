module GTOC12
    using DataFrames

    export asteroid_df, planet_df, ET₀

    # File paths
    const PROBLEM_DATA = joinpath(@__DIR__, "../problem")

    # Dynamics
    include("dynamics/keplerian_elements.jl")
    include("dynamics/modified_equinoctial_elements.jl")
    include("dynamics/universal_variable.jl")
    include("dynamics/frame_conversions.jl")
    include("dynamics/naive_flyby.jl")


    # Utilities
    include("utils/SPICE_wrappers.jl")
    include("utils/constants.jl")
    include("utils/asteroids.jl")
    include("utils/defaults.jl")
    include("utils/helper_functions.jl")
    include("utils/planets_on_rails.jl")

    # Targeters and Optimizers
    include("targeting/targeting.jl")

    # Plotting
    include("plotting/plotting.jl")

    function __init__()
        furnish_all_kernels()
        # Only load/manipulate dataframes once
        global asteroid_df = get_asteroid_df()
        global planet_df = get_planet_df()

        # Problem start time
        global ET₀ = asteroid_df[1, :ET]
    end

end
