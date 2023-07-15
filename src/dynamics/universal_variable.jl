export propagate_universal

# Vallado Algorithm 1
function get_c2_c3(ψ)
    if ψ > 1e-6
        c2 = (1 - cos(√ψ)) / ψ
        c3 = (√ψ - sin(√ψ)) / √(ψ^3)
    else
        if ψ < -1e-6
            c2 = (1 - cosh(√(-ψ))) / ψ
            c3 = (sinh(√(-ψ)) - √(-ψ)) / √(-ψ^3)
        else
            c2 = 1.0 / 2.0
            c3 = 1.0 / 6.0
        end
    end
    return c2, c3
end

# Vallado Algorithm 8
function solve_for_universal_variable(; Χ₀, α, μ=GTOC12.μ_☉, r⃗₀, v⃗₀, Δt, tol=1e-8, MAX_ITER=10 )
    r₀ = norm(r⃗₀)

    Χₙ = Χ₀
    num_iter = 0
    while true
        ψ = Χₙ^2 * α

        c2, c3 = get_c2_c3(ψ)

        r = Χₙ^2 * c2 + (r⃗₀ ⋅ v⃗₀) / √μ * Χₙ * (1 - ψ * c3) + r₀ * (1 - ψ * c2)

        Χₙ⁺ = Χₙ + (√μ * Δt - Χₙ^3 * c3 - (r⃗₀ ⋅ v⃗₀) / √μ * Χₙ^2 * c2 - r₀ * Χₙ * (1 - ψ * c3)) / r

        # Check if we meet break conditions
        abs(Χₙ⁺ - Χₙ) > tol || return Χₙ⁺, ψ, r

        # Add error for violation of max iterations
        if num_iter > MAX_ITER
            err = ErrorException("Universal Variable Propagation: Max number of iterations reached in Kepler's Equation. Check your units.")
            throw(err)
        end

        Χₙ = Χₙ⁺

    end

end

# Vallado Algorithm 8
function propagate_universal(x⃗, Δt; μ=GTOC12.μ_☉, tol=1e-8)
    r⃗₀ = x⃗[1:3]
    v⃗₀ = x⃗[4:6]

    # Need to normalize by canonical units
    CDU = norm(r⃗₀)  # Canonical Distance Unit
    CTU = √(CDU^3/μ) # Canonical Time Unit
    μ = 1

    r⃗₀ = r⃗₀/CDU
    v⃗₀ = v⃗₀/(CDU/CTU)
    Δt = Δt/CTU


    r₀ = norm(r⃗₀)
    v₀ = norm(v⃗₀)
    ξ = v₀^2 / 2.0 - μ / r₀
    α = -v₀^2 / μ + 2.0 / r₀

    # Initial guess at universal variable based on regime
    if α > 0.000001 # Ellipse
        Χ₀ = √μ * Δt * α
    elseif abs(α) < 0.000001 # Parabola
        h⃗ = r⃗₀ × v⃗₀
        p = norm(h⃗)^2 / μ
        s = acot(3 / 2 * √(μ / p^3) * Δt)
        w = atan(tan(s)^(1/3))
        Χ₀ = √p * 2 * cot(2w)
    else # hyperbola
        a = 1 / α
        Χ₀ = sign(Δt) * √(-a) * log(-2μ * α * Δt / (r⃗₀ ⋅ v⃗₀ + sign(Δt) * √(-μ*a) * 1 - r₀ * α))
    end

    # Solve for the universal variable
    Χₙ⁺, ψ, r = solve_for_universal_variable(; Χ₀=Χ₀, Δt=Δt, α=α, μ=μ, r⃗₀=r⃗₀, v⃗₀=v⃗₀, tol=tol )

    # Calculate f, g, ḟ, ġ, and resultant Cartesian state
    c2, c3 = get_c2_c3(ψ)

    f = 1 - Χₙ⁺^2 / r₀ * c2
    g = Δt - Χₙ⁺^3 / √μ * c3
    ġ = 1 - Χₙ⁺^2 / r * c2
    ḟ = √μ / (r * r₀) * Χₙ⁺ * (ψ * c3 - 1)

    r⃗ = f * r⃗₀ + g * v⃗₀
    v⃗ = ḟ * r⃗₀ + ġ * v⃗₀

    return vcat(r⃗*CDU, v⃗*CDU/CTU) # Re=normalize and return
end



