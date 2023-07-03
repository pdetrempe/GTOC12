module GTOC12
    using DataFrames
    using CSV
    using SPICE
    using StaticArrays
    using LinearAlgebra
    using ForwardDiff, FiniteDiff
    using RobotDynamics
    using Altro
    using TrajectoryOptimization


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
    include("utils/celestial_body.jl")
    include("utils/asteroids.jl")
    include("utils/planets.jl")
    include("utils/defaults.jl")
    include("utils/helper_functions.jl")

    # Targeters and Optimizers
    include("targeting/targeting.jl")
    include("targeting/lambert_solver.jl")
    include("targeting/spacecraft.jl")
    include("optimization/optimization.jl")

    # Plotting
    include("plotting/plotting.jl")








end
