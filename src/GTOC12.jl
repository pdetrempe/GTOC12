module GTOC12
    using DataFrames


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

    # Targeters and Optimizers
    include("targeting/targeting.jl")
    include("targeting/lambert_solver.jl")

    # Plotting
    include("plotting/plotting.jl")






end
