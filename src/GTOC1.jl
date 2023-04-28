module GTOC1
    using Reexport

    @reexport using SPICE
    @reexport using AstroTime
    @reexport using Plots
    @reexport using PlanetOrbits
    @reexport import PlanetOrbits: m2au, _trueanom_from_eccanom
    @reexport using OrdinaryDiffEq
    @reexport using LinearAlgebra

    # Dynamics
    include("dynamics/keplerian_elements.jl")
    include("dynamics/modified_equinoctial_elements.jl")

    # Utilities
    include("utils/time.jl")
    include("utils/SPICE_kernels.jl")

end
