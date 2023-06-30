using LinearAlgebra

export lambert

# Final bit of Vallado Algorithm 58
function get_lambert_velocities(; r⃗₀, r⃗, A, yₙ, μ)
    r₀ = norm(r⃗₀)
    r = norm(r⃗)

    f = 1 - yₙ / r₀
    ġ = 1 - yₙ / r
    g = A * √(yₙ / μ)

    v⃗₀ = (r⃗ - f * r⃗₀) / g
    v⃗ = (ġ * r⃗ - r⃗₀) / g

    return v⃗₀, v⃗
end

# Vallado Algorithm 58
# NOTE: This algorithm is intended to work with Canonical units (notice the default μ value)
# Greatly simplifying the equations and increasing numerical stability
#
# tₘ is +1 for short way transfers and -1 for long way.
# TODO: Add Enum for tₘ
function lambert(; r⃗₀, r⃗, Δt, tₘ, μ=GTOC12.μ_☉, MAX_ITER=50, tol=1e-6)

    # Normalize by canonical units
    CDU = norm(r⃗₀)  # Canonical Distance Unit
    CTU = √(CDU^3 / μ)# Canonical Time Unit
    μ = 1

    r⃗₀ = r⃗₀ ./ CDU
    r⃗ = r⃗ ./ CDU
    Δt = Δt ./ CTU

    r₀ = norm(r⃗₀)
    r = norm(r⃗)
    cosΔν = r⃗₀ ⋅ r⃗ / (r₀ * r)
    sinΔν = tₘ * √(r * r₀ * (1 - cosΔν^2))
    A = tₘ * √(r * r₀ * (1 + cosΔν))

    if abs(A) < eps(Float64)
        err = "Invalid Lambert arc angle"
        throw(err)
    end

    ψₙ = 0.0
    c2 = 1 / 2
    c3 = 1 / 6

    # Loop until time convergence
    ψ_up, ψ_low = 4 * π^2, -4 * π^2
    num_iter = 0
    while true
        yₙ = r₀ + r + A * (ψₙ * c3 - 1) / √c2

        if A > 0 && yₙ < 0
            # readjust ψ_low until y > 0
            MUST_ADJUST_ψ_low = true
        else
            MUST_ADJUST_ψ_low = false
        end

        Χₙ = √(yₙ / c2)
        Δtₙ = (Χₙ^3 * c3 + A * √yₙ) / √μ

        if Δtₙ <= Δt || MUST_ADJUST_ψ_low
            ψ_low = ψₙ
        else
            ψ_up = ψₙ
        end

        ψₙ⁺ = (ψ_up + ψ_low) / 2
        println("ψₙ⁺ = "*string(ψₙ⁺))

        c2, c3 = get_c2_c3(ψₙ⁺)

        ψₙ = ψₙ⁺

        # Check if we meet break conditions
        # If so, calculate velocities, re-dimensionalize, and return
        abs(Δtₙ - Δt) > tol || return get_lambert_velocities(; r⃗₀=r⃗₀, r⃗=r⃗, A=A, yₙ=yₙ, μ=μ) .* CDU ./ CTU

        # Add error for violation of max iterations
        if num_iter > MAX_ITER
            err = ErrorException("lambert_solver: Max number of iterations reached in Lambert Solver. Check your units.")
            throw(err)
        end
        num_iter += 1
        println(num_iter)
    end
end