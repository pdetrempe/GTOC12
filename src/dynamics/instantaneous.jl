
## Launch vehicle impulse
LV_condition(x, t, integrator) = t == 0

function LV_affect!(integrator)
    MEE = integrator.u        # State is Modified Equinoctal Elements
    μ = integrator.p.μ
    ΔV = integrator.p.ΔV_LV_inrt # Downside: DV has to be predefined

    x⃗₋ = MEE2Cartesian(MEE; μ=μ)

    x⃗₊ = x⃗₋ + ΔV 

    integrator.u = Cartesian2MEE(x⃗₊; μ=μ) # State is Modified Equinoctal Elements

end

LV_callback = DiscreteCallback(LV_condition, LV_affect!)