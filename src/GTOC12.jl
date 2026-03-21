module GTOC12
using DataFrames
using CSV
using SPICE
using StaticArrays
using LinearAlgebra
using ForwardDiff, FiniteDiff
# using RobotDynamics  # Commented out - breaks Julia 1.12
# using Altro  # Commented out - depends on RobotDynamics
# using TrajectoryOptimization  # Commented out - depends on RobotDynamics
# using OrdinaryDiffEq  # Commented out - pulls in DifferentialEquations which breaks

# File paths
const PROBLEM_DATA = joinpath(@__DIR__, "../problem")

# Dynamics
include("dynamics/keplerian_elements.jl")
include("dynamics/modified_equinoctial_elements.jl")
include("dynamics/universal_variable.jl")
include("dynamics/frame_conversions.jl")


# Utilities
include("utils/SPICE_wrappers.jl")
include("utils/constants.jl")
include("utils/celestial_body.jl")
include("utils/asteroids.jl")
include("utils/planets.jl")
include("dynamics/naive_flyby.jl")
include("utils/defaults.jl")
include("utils/helper_functions.jl")
include("utils/structures.jl")
include("utils/mining_functions.jl")
include("utils/data_recording_functions.jl")
include("utils/planets_on_rails.jl")

# Targeters, Optimizers, and Controllers
include("targeting/targeting.jl")
include("targeting/lambert_solver.jl")
# include("targeting/spacecraft.jl")  # Disabled - requires RobotDynamics
# include("optimization/optimization.jl")  # Disabled - requires Altro/TrajectoryOptimization
include("optimization/infiniteopt_low_thrust.jl")  # This one works standalone!
include("dynamics/primer_vector.jl")
# include("controls/controls.jl")  # Disabled - DifferentialEquations broken

# Plotting
include("plotting/plotting.jl")


end
