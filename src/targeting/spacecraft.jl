using StaticArrays
using RobotDynamics: @autodiff
const RD = RobotDynamics

export ImpulsiveSpacecraft, RD

@autodiff struct ImpulsiveSpacecraft <: RD.DiscreteDynamics
    ET_launch::Float64
    target_asteroid::Asteroid
    μ::Float64
end

function ImpulsiveSpacecraft(;ET_launch::Float64, target_asteroid::Asteroid, μ::Float64=GTOC12.μ_☉)
    return ImpulsiveSpacecraft(ET_launch, target_asteroid, μ)
end

function Base.copy(c::ImpulsiveSpacecraft)
    ImpulsiveSpacecraft(c.ET_launch, c.target_asteroid, c.μ)
end

# In-place dynamics
function RD.discrete_dynamics!(model::ImpulsiveSpacecraft, ẋ, x, u, t, dt)
    ΔV = u
    x₀⁺ = x + [0;0;0; ΔV[:] ]
    ẋ[:] =  propagate_universal(x₀⁺, dt; μ=model.μ)
    nothing
end

# Out-of-place dynamics
function RD.discrete_dynamics(model::ImpulsiveSpacecraft, x, u, t, dt)
    ΔV = u
    x₀⁺ = x + [0;0;0; ΔV[:] ]
    xₜ =  propagate_universal(x₀⁺, dt; μ=model.μ)
    SVector{length(xₜ)}(xₜ)
end

RD.state_dim(::ImpulsiveSpacecraft) = 6
RD.control_dim(::ImpulsiveSpacecraft) = 3