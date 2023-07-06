
export optimize_impulsive_transfer, get_Φ, get_lambert_ΔV_and_p⃗

# 1. Use autodiff of propagate_universal to get STM
function get_Φ(x₀, Δt; μ=GTOC12.μ_☉)
    Φ = ForwardDiff.jacobian( x -> propagate_universal(x, Δt; μ), x₀ )
end

# This follows chapter 5 in Prussing-Optimal Spacecraft Trajectories
get_p₀(ΔV₀) = ΔV₀/norm(ΔV₀) # Eq 5.3
get_p_f(ΔVₜ) = ΔVₜ/norm(ΔVₜ)

function propagate_p_and_ṗ(p₀, ṗ₀, Φ) # Eq 5.5
    # Used to propagate primer vectore and derivative
    xₜ = Φ * [p₀, ṗ₀]
    pₜ = xₜ[1:3]
    ṗₜ = xₜ[4:6]
    return pₜ, ṗₜ
end

M(Φ) = Φ[1:3, 1:3]
N(Φ) = Φ[1:3, 4:6]
S(Φ) = Φ[4:6, 1:3]
T(Φ) = Φ[4:6, 4:6]

# Prussing Eq. 5.10
function get_primer_vector(p₀, p_f, Φ_0_to_t, Φ_0_to_f)
    N_t0 = N(Φ_0_to_t)
    N_f0 = N(Φ_0_to_f)
    M_t0 = M(Φ_0_to_t)
    M_f0 = M(Φ_0_to_f)

    pₜ = N_t0/N_f0\p_f + [M_t0 - N_t0/N_f0\M_f0]*p₀
end

# Prussing 5.9
function get_ṗ₀(p₀, p_f, Φ_0_to_f)
    N_f0 = N(Φ_0_to_f)
    M_f0 = M(Φ_0_to_f)

    return N_f0\(p_f - M_f0*p₀)
end

# Prussing 5.8
function get_ṗ_f(p₀, ṗ₀, Φ_0_to_f)
    S_f0 = S(Φ_0_to_f)
    T_f0 = T(Φ_0_to_f)

    return S_f0*p₀ + T_f0*ṗ₀
end

# Function to get cost from bodies and Lambert arc
function get_lambert_ΔV_and_p⃗(t0, tf; body_from::CelestialBody, body_to::CelestialBody, ET_start)
    Δt = tf - t0

    # Get initial/final states
    x₀ = get_body_state(body_from; ET=ET_start+t0)
    ET_target = ET_start + tf
    x_target = get_body_state(body_to; ET=ET_target)

    # # Non-dimensionalize everything
    # x₀, CDU, CTU, μ_canonical = get_canonical_state(x₀; μ=GTOC12.μ_☉)
    # x_target = get_canonical_state(x_target, CDU, CTU)
    # Δt = canonical_time(Δt; CTU=CTU)


    # Find Lambert transfer between two positions
    # TODO: Make a lambert type?
    println("Δt = $Δt")
    v⃗₀⁺, v⃗ = lambert(; r⃗₀=x₀[1:3], r⃗=x_target[1:3], Δt=Δt, tₘ=1) # μ=μ_canonical)
    ΔV₀ = v⃗₀⁺ - x₀[4:6]
    ΔVₜ = x_target[4:6] - v⃗
    x₀⁺ = x₀ + [0;0;0;v⃗₀⁺[:]]
    println("x₀ = $x₀")
    println("x₀⁺ = $x₀⁺")

    # Get STM from Lambert arc
    Φ = get_Φ(x₀⁺, Δt) #; μ=μ_canonical)
        
    # Get Primer vector conditions
    p₀ = get_p₀(ΔV₀) # Eq 5.3
    p_f = get_p_f(ΔVₜ)
    ṗ₀ = get_ṗ₀(p₀, p_f, Φ)
    ṗ_f = get_ṗ_f(p₀, ṗ₀, Φ)

    # # Redimensionalize everything
    # redimensionalize_state!(x₀; CDU=CDU, CTU=CTU)
    # redimensionalize_state!(x_target; CDU=CDU, CTU=CTU)
    # # redimensionalize_vel!(ΔV₀; CDU=CDU, CTU=CTU)
    # # redimensionalize_vel!(ΔVₜ; CDU=CDU, CTU=CTU)


    return ΔV₀, ΔVₜ, p₀, ṗ₀, p_f, ṗ_f
end

# Prussing 5.37
function dJ_dtf(t0, tf; body_from::CelestialBody, body_to::CelestialBody, ET_start)
    ΔV₀, ΔVₜ, p₀, ṗ₀, p_f, ṗ_f = get_lambert_ΔV_and_p⃗(t0, tf; body_from=body_from, body_to=body_to, ET_start=ET_start)

    -norm(ΔVₜ)*ṗ_f'*p_f
end


# Prussing 5.36
function dJ_dt0(t0, tf; body_from::CelestialBody, body_to::CelestialBody, ET_start)
    ΔV₀, ΔVₜ, p₀, ṗ₀, p_f, ṗ_f = get_lambert_ΔV_and_p⃗(t0, tf; body_from=body_from, body_to=body_to, ET_start=ET_start)

    -norm(ΔV₀)*ṗ₀'*p₀
end

function newton_step(x, f, ḟ; newton_step=1)
    x⁺ = x - newton_step * f/ḟ
end

# A few approaches here:
# Given explicit function for J = ΣΔV ( Lambert(t0, tf, Δt) ), can try just finding the Hessian
# Or use analytical functions for dJ from Prussing primer vector and:
#   - Do a line search to find where conditions met. Iterate on gradient w.r.t. start/final times
#   - Take Jacobian of dJ function to get ddJ for Newton's method (I think this is cooler)

scalar_ṗ(p⃗, p⃗̇) = p⃗̇'*p⃗

function optimize_impulsive_transfer(body_from::CelestialBody, body_to::CelestialBody; ET_start, Δt_guess, bound_initial_time=false, bound_final_time=false, tol=1e-8, MAX_ITER=20)
    t0 = 0
    tf = Δt_guess
    ΔV₀, ΔVₜ, p⃗₀, p⃗̇₀, p⃗_f, p⃗̇_f= get_lambert_ΔV_and_p⃗(t0, tf; body_from=body_from, body_to=body_to, ET_start=ET_start)

    # TODO: Make this conditional check EITHER initial/final time, depending on boolean inputs above
    # Can maybe pass a function to a variable and use that in the while conditional
    ṗ₀ = scalar_ṗ(p⃗₀, p⃗̇₀)
    ṗ_f = scalar_ṗ(p⃗_f, p⃗̇_f)

    println("p⃗₀ = $p⃗₀")
    println("p⃗̇₀ = $p⃗̇₀")
    println("p⃗_f = $p⃗_f")
    println("p⃗̇_f = $p⃗̇_f")


    num_iter = 0
    while abs(ṗ₀) > tol && abs(ṗ_f) > tol
        J = ΔV₀ + ΔVₜ # total cost
        print("J = $J")
        ∂J_∂t0 = dJ_dt0(t0, tf; body_from=body_from, body_to=body_to, ET_start=ET_start)
        ∂J_∂tf = dJ_dtf(t0, tf; body_from=body_from, body_to=body_to, ET_start=ET_start)


        # Add a function call for DJ_dt0 as a function of all the above
        ∂2J_dt₀2 = ForwardDiff.derivative( t0 -> dJ_dt0(t0, tf; body_from, body_to, ET_start ), t0 )
        ∂2J_dtf2 = ForwardDiff.derivative( tf -> dJ_dtf(t0, tf; body_from, body_to, ET_start ), tf )

        # ForwardDiff ^ That function call to get the 2nd derivative w.r.t. t0/tf

        # Take a Newton step and update t0 (ET_start + dt0), tf (ET_start + Δt + Δtf)
        println("t0 = $t0")
        println("tf = $tf")
        partial_newton = 0.1

        t0 = newton_step(t0, ∂J_∂t0, ∂2J_dt₀2; newton_step = partial_newton) 
        tf = newton_step(tf, ∂J_∂tf, ∂2J_dtf2; newton_step = partial_newton) 
        println("t0⁺ = $t0")
        println("tf⁺ = $tf")

        # Resolve for DVs and primer vector values
        ΔV₀, ΔVₜ, p⃗₀, p⃗̇₀, p⃗_f, p⃗̇_f = get_lambert_ΔV_and_p⃗(t0, tf; body_from=body_from, body_to=body_to, ET_start=ET_start)
        ṗ₀ = scalar_ṗ(p⃗₀, p⃗̇₀)
        ṗ_f = scalar_ṗ(p⃗_f, p⃗̇_f)

        # Add error for violation of max iterations
        println("$num_iter")
        num_iter = num_iter + 1
        if num_iter > MAX_ITER
            err = ErrorException("Primer vector Newton Method: Max number of iterations reached in Primer vector Newton Method. Check your units.")
            throw(err)
        end
    end


    return ΔV₀, ΔVₜ, t0, tf
end
