using StaticArrays
using RobotDynamics: @autodiff
const RD = RobotDynamics

export ImpulsiveSpacecraft, RD, ContinuousSpacecraft, Continuous2Spacecraft

#------------------------------- Impulsive Spacecraft
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




# ------------------------------- Continuous Spacecraft ------------------------------- #
@autodiff struct ContinuousSpacecraft <: RD.ContinuousDynamics
    ET_launch::Float64
    target_asteroid::Asteroid
    μ::Float64
    #T::Float64 
    #Isp::Float64
    # TODO: Add the thrust/Isp/etc. values in here
end

function ContinuousSpacecraft(;ET_launch::Float64, target_asteroid::Asteroid, μ::Float64=GTOC12.μ_☉)
    return ContinuousSpacecraft(ET_launch, target_asteroid, μ)
end

function Base.copy(c::ContinuousSpacecraft)
    ContinuousSpacecraft(c.ET_launch, c.target_asteroid, c.μ)
end

# Add continuous dynamics functions
function calculate_hamiltonian(x, p) 
    u, μ = p

    MEE = x[1:6]
    m = x[7]
    λ = x[8:13]

    A = A_equinoctial(MEE; μ=μ)
    B = B_equinoctial(MEE; μ=μ)

    # The below works since norm(u) = T*δ_star/m
    H = norm(u/m) +    λ'*A    +   λ'*B*u
    return H
end

#function get_opt_control_from_state(model::ContinuousSpacecraft, state)
#    MEE = state[1:6]; m = state[7];  #λ = state[8:13]
#
#    # Solve entire problem in MEE
#    # B = B_equinoctial(MEE; μ=model.μ)
#
#    # # Use this function during "roll-out" for initial guess of control values
#    # u_star = -B'*λ/norm(B'*λ)  
#    # S = norm(B'*λ) .- 1        
#    # if S > 0
#    #     δ_star = 1
#    # else
#    #     δ_star = 0
#    # end
#    # u = T_max*δ_star*u_star
#    u = [0.,0.,0.]
#end

function RD.dynamics(model::ContinuousSpacecraft, x, u, t)
    MEE = x[1:6]; m = x[7]; # λ = x[8:13]

    # Solve entire problem in MEE
    A = A_equinoctial(MEE; μ=model.μ)
    B = B_equinoctial(MEE; μ=model.μ)

    # Spacecraft Dynamics 
    # *************************************************
    # u = T*δ_star/m * u_star
    dx = A + B*u/m

    if norm(u) > 0
        δ_star=1
    else
        δ_star=0
    end

    # Mass variation
    # *************************************************
    # exhaust velocity
    c = GTOC12.Isp * GTOC12.g0 # specific impulse and grav accel at sea level (m/s)
    dm = -GTOC12.T_max/c * δ_star

    # Costate diff eqs
   # p = u, model.μ
   # dH_dx = ForwardDiff.gradient(x -> calculate_hamiltonian(x, p ), x) #[1:6])
   # dλ = -dH_dx'

    dstate = [dx[:]; dm] #; dλ[1:6]]
    SVector{length(dstate)}(dstate)

end

function RD.dynamics!(model::ContinuousSpacecraft, ẋ, x, u, t)
    ẋ[:] = RD.dynamics(model, x, u, t)
    nothing
end



RD.state_dim(::ContinuousSpacecraft) = 7 #13
RD.control_dim(::ContinuousSpacecraft) = 3



# ------------------------------- Continuous2 Spacecraft ------------------------------- #
#@autodiff struct Continuous2Spacecraft <: RD.ContinuousDynamics
#    ET_launch::Float64
#    target_asteroid::Asteroid
#    μ::Float64
#    #T::Float64 
#    #Isp::Float64
#    # TODO: Add the thrust/Isp/etc. values in here
#end
#
#function Continuous2Spacecraft(;ET_launch::Float64, target_asteroid::Asteroid, μ::Float64=GTOC12.μ_☉)
#    return Continuous2Spacecraft(ET_launch, target_asteroid, μ)
#end
#
#function Base.copy(c::Continuous2Spacecraft)
#    Continuous2Spacecraft(c.ET_launch, c.target_asteroid, c.μ)
#end
#
## Add continuous dynamics functions
#function calculate_hamiltonian(x, p) 
#    u, μ = p
#
#    MEE = x[1:6]
#    m = x[7]
#    λ = x[8:13]
#
#    A = A_equinoctial(MEE; μ=μ)
#    B = B_equinoctial(MEE; μ=μ)
#
#    # The below works since norm(u) = T*δ_star/m
#    H = norm(u/m) +    λ'*A    +   λ'*B*u
#    return H
#end
#
#function get_opt_control_from_state(model::Continuous2Spacecraft, state)
#    MEE = state[1:6]; m = state[7]; λ = state[8:13]
#
#    # Solve entire problem in MEE
#    B = B_equinoctial(MEE; μ=model.μ)
#
#    # Use this function during "roll-out" for initial guess of control values
#    u_star = -B'*λ/norm(B'*λ)  
#    S = norm(B'*λ) .- 1        
#    if S > 0
#        δ_star = 1
#    else
#        δ_star = 0
#    end
#    u = T_max*δ_star*u_star
#end
#
#function RD.dynamics(model::Continuous2Spacecraft, x, u, t)
#    MEE = x[1:6]; m = x[7]; λ = x[8:13]
#
#    # Solve entire problem in MEE
#    A = A_equinoctial(MEE; μ=model.μ)
#    B = B_equinoctial(MEE; μ=model.μ)
#
#    # Spacecraft Dynamics 
#    # *************************************************
#    # u = T*δ_star/m * u_star
#    dx = A + B*u/m
#
#    if norm(u) > 0
#        δ_star=1
#    else
#        δ_star=0
#    end
#
#    # Mass variation
#    # *************************************************
#    # exhaust velocity
#    c = GTOC12.Isp * GTOC12.g0 # specific impulse and grav accel at sea level (m/s)
#    dm = -GTOC12.T_max/c * δ_star
#
#    # Costate diff eqs
#    p = u, model.μ
#    dH_dx = ForwardDiff.gradient(x -> calculate_hamiltonian(x, p ), x) #[1:6])
#    dλ = -dH_dx'
#
#    dstate = [dx[:]; dm; dλ[1:6]]
#    SVector{length(dstate)}(dstate)
#
#end
#
#function RD.dynamics!(model::Continuous2Spacecraft, ẋ, x, u, t)
#    ẋ[:] = RD.dynamics(model, x, u, t)
#    nothing
#end
#
#
#
#RD.state_dim(::Continuous2Spacecraft) = 13
#RD.control_dim(::Continuous2Spacecraft) = 3
