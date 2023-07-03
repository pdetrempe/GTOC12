using GTOC12
using Plots
using LinearAlgebra

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1))
asteroid = GTOC12.asteroid_df[ID_min, :]

# Try out shooting method to hit asteroid
furnish_all_kernels()
DV₀ = [-1000; 0; 0]
x₀ = get_planet_state("EARTH", GTOC12.ET₀) + [0; 0; 0; DV₀[:]]

# Transfer time
Δt = 1 / 2 * 365 * 24 * 3600
ΔV∞_max = 6000.0;

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = GTOC12.ET₀ + Δt
x_target = get_asteroid_state(asteroid, ET_target)
r_target = view(x_target, 1:3)

# NOTE: This is super brittle to initial conditions
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
x₀⁺, xₜ = fixed_time_single_shoot(x₀, Δt, r_target; print_iter=false)
ΔV₀ = (x₀⁺-x₀)[4:6]
ΔVₜ = (xₜ-x_target)[4:6]

#--------------------- Plot Earth/asteroid/spacecraft
# ETs = [GTOC12.ET₀, ET_target]
# plot_asteroid_from_df_row(asteroid; ET_in=ETs)
# plot_planet!(planet="EARTH"; ET_in=ETs, label="Earth", color=colormap("Blues"))
# plot_coast!(x₀⁺, Δt; label="Coast 1", color=colormap("Greens"))

# 1. Use autodiff of propagate_universal to get STM

# This follows chapter 5 in Prussing-Optimal Spacecraft Trajectories
p₀ = ΔV₀/norm(ΔV₀) # Eq 5.3
pₜ = ΔVₜ/norm(ΔVₜ)

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

    pₜ = N_t0*inv(N_f0)*p_f + [M_t0 - N_t0*inv(N_f0)*M_f0]*p₀
end

# Prussing 5.9
function get_ṗ₀(p₀, p_f, Φ_0_to_f)
    N_f0 = N(Φ_0_to_f)
    M_f0 = M(Φ_0_to_f)

    return inv(N_f0)*(p_f - M_f0*p₀)
end

# Prussing 5.8
function get_ṗ_f(p₀, ṗ₀, Φ_0_to_f)
    S_f0 = S(Φ_0_to_f)
    T_f0 = T(Φ_0_to_f)

    return S_f0*p₀ + T_f0*ṗ₀
end

# function dJ(;ΔV₀, ΔVₜ, p₀, ṗ₀, pₜ, ṗₜ, dt₀, dtₜ)

# end

# A few approaches here:
# Given explicit function for J = ΣΔV ( Lambert(t0, tf, Δt) ), can try just finding the Hessian
# Or use analytical functions for dJ from Prussing primer vector and:
#   - Do a line search to find where conditions met. Iterate on gradient w.r.t. start/final times
#   - Take Jacobian of dJ function to get ddJ for Newton's method (I think this is cooler)

function optimize_transfer_times(p₀, p_f, Φ_0_to_f; bound_initial_time=false, bound_final_time=false, tol=1e-6)
    ṗ₀ = get_ṗ₀(p₀, p_f, Φ_0_to_f)
    ṗ_f = get_ṗ_f(p₀, ṗ₀, Φ_0_to_f)

    # TODO: Make this conditional check EITHER initial/final time, depending on boolean inputs above
    # Can maybe pass a function to a variable and use that in the while conditional

    while abs(ṗ₀) > tol && abs(ṗ_f) > tol
        
    end
    

    # See page 49 of Prussing for these conditionals
    # if ṗ₀ > 0
    #     ADD_INITIAL_COAST = true
    # end

    # if ṗ_f < 0
    #     ADD_FINAL_COAST = true
    # end


end
