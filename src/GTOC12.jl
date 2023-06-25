module GTOC12
    using SPICE
    using AstroTime
    # using Plots
    # using PlanetOrbits
    # import PlanetOrbits: m2au, _trueanom_from_eccanom
    # using OrdinaryDiffEq
    # using LinearAlgebra


    # File paths
    const PROBLEM_DATA = joinpath(@__DIR__, "../problem")

    # Dynamics
    include("dynamics/keplerian_elements.jl")
    include("dynamics/modified_equinoctial_elements.jl")
    include("dynamics/frame_conversions.jl")
    include("dynamics/instantaneous.jl")
    include("dynamics/naive_flyby.jl")


    # Utilities
    include("utils/SPICE_kernels.jl")
    include("utils/asteroids.jl")
    include("utils/constants.jl")
    include("utils/time.jl")
    include("utils/defaults.jl")




end
