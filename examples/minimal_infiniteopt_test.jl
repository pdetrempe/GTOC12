"""
Minimal InfiniteOpt Low-Thrust Test

Tests the InfiniteOpt module with minimal dependencies.
This bypasses the full GTOC12 module to avoid compatibility issues.
"""

# Minimal imports
using InfiniteOpt
using Ipopt
using LinearAlgebra
using Statistics
using Printf

println("="^70)
println("Minimal InfiniteOpt Low-Thrust Test")
println("="^70)

# Constants (from GTOC12)
const μ_☉ = 1.32712440018e11 * (1000)^3  # m³/s²
const AU_m = 1.49597870691e11  # m
const T_max = 0.6  # N
const Isp = 4000.0  # s
const g0 = 9.80665  # m/s²

# Simple test problem: circular orbit transfer
println("\nTest Problem: Simple Orbit Transfer")
println("-"^70)

# Initial state: Circular orbit at 1 AU
r0_mag = 1.0 * AU_m
v0_mag = sqrt(μ_☉ / r0_mag)
x0 = [r0_mag; 0.0; 0.0; 0.0; v0_mag; 0.0; 2200.0]  # [r; v; m] in m, m/s, kg
println("Initial orbit radius: $(r0_mag / AU_m) AU")
println("Initial velocity: $(v0_mag / 1000) km/s")
println("Initial mass: $(x0[7]) kg")

# Target state: Circular orbit at 1.1 AU
rf_mag = 1.1 * AU_m
vf_mag = sqrt(μ_☉ / rf_mag)
xf = [rf_mag; 0.0; 0.0; 0.0; vf_mag; 0.0]  # [r; v] in m, m/s
println("Target orbit radius: $(rf_mag / AU_m) AU")
println("Target velocity: $(vf_mag / 1000) km/s")

# Transfer time (increased for feasibility)
t_f = 200.0 * 24 * 3600  # 200 days in seconds
println("Transfer time: $(t_f / (24*3600)) days")

# Discretization
N = 50  # Increased for better accuracy
println("Discretization points: $N")

println("\n" * "="^70)
println("Running InfiniteOpt Optimization")
println("="^70)

# =============================================================================
# CANONICAL UNITS
# =============================================================================
L_scale = AU_m
T_scale = 10.0 * 86400.0
V_scale = L_scale / T_scale
M_scale = x0[7]
F_scale = M_scale * L_scale / T_scale^2

# Non-dimensionalize
x0_nd = zeros(7)
x0_nd[1:3] = x0[1:3] / L_scale
x0_nd[4:6] = x0[4:6] / V_scale
x0_nd[7] = x0[7] / M_scale

xf_nd = zeros(6)
xf_nd[1:3] = xf[1:3] / L_scale
xf_nd[4:6] = xf[4:6] / V_scale

t_f_nd = t_f / T_scale
μ_nd = μ_☉ / (L_scale^3 / T_scale^2)
T_max_nd = T_max / F_scale
c_dim = Isp * g0
c_nd = c_dim / V_scale
m_min_nd = 500.0 / M_scale

println("Canonical scales:")
println("  L_scale = $(L_scale / AU_m) AU")
println("  T_scale = $(T_scale / 86400) days")
println("  V_scale = $(V_scale / 1000) km/s")
println("  M_scale = $(M_scale) kg")

# =============================================================================
# BUILD INFINITEOPT MODEL
# =============================================================================
println("\nBuilding InfiniteOpt model...")

model = InfiniteModel(optimizer_with_attributes(
    Ipopt.Optimizer,
    "print_level" => 5,
    "max_iter" => 3000,
    "tol" => 1e-5,
    "constr_viol_tol" => 1e-4,
    "acceptable_tol" => 1e-3,
    "acceptable_iter" => 15,
    "mu_strategy" => "adaptive",
    "max_cpu_time" => 600.0
))

# Time parameter
@infinite_parameter(model, t in [0, t_f_nd], num_supports = N,
    derivative_method = InfiniteOpt.OrthogonalCollocation(4))

# State variables
@variable(model, r[1:3], Infinite(t))
@variable(model, v[1:3], Infinite(t))
@variable(model, m, Infinite(t))

# Control variables (spherical angles + throttle)
@variable(model, 0 <= θ <= π, Infinite(t))
@variable(model, φ, Infinite(t))
@variable(model, 0 <= δ <= 1, Infinite(t))

# Piecewise constant controls (optional - may not be available in this version)
# Try to make controls piecewise constant if function is available
try
    InfiniteOpt.constant_over_collocation(θ, t)
    InfiniteOpt.constant_over_collocation(φ, t)
    InfiniteOpt.constant_over_collocation(δ, t)
    println("✓ Using piecewise constant controls")
catch e
    println("⚠️  Piecewise constant controls not available - using continuous controls")
    println("   (This is OK for testing, just less efficient)")
end

# Initial guesses (linear interpolation)
for i in 1:3
    set_start_value_function(r[i], t_val -> x0_nd[i] + (xf_nd[i] - x0_nd[i]) * t_val / t_f_nd)
    set_start_value_function(v[i], t_val -> x0_nd[i+3] + (xf_nd[i+3] - x0_nd[i+3]) * t_val / t_f_nd)
end
set_start_value_function(m, t_val -> x0_nd[7] * (1.0 - 0.1 * t_val / t_f_nd))
set_start_value_function(θ, t_val -> π / 2)
set_start_value_function(φ, t_val -> 0.0)
set_start_value_function(δ, t_val -> 0.5)

# Objective: Maximize final mass
t_final = upper_bound(t)
@objective(model, Max, m(t_final))

# Initial conditions
@constraint(model, [i = 1:3], r[i](0) == x0_nd[i])
@constraint(model, [i = 1:3], v[i](0) == x0_nd[i+3])
@constraint(model, m(0) == x0_nd[7])

# Final conditions (intercept only for this test)
@constraint(model, [i = 1:3], r[i](t_final) == xf_nd[i])

# Dynamics: dr/dt = v
@constraint(model, [i = 1:3], ∂(r[i], t) == v[i])

# Dynamics: dv/dt = a_grav + a_thrust
for i in 1:3
    if i == 1
        @constraint(model,
            ∂(v[1], t) == -μ_nd * r[1] / (r[1]^2 + r[2]^2 + r[3]^2)^(3 / 2)
                          + (T_max_nd / m) * δ * sin(θ) * cos(φ)
        )
    elseif i == 2
        @constraint(model,
            ∂(v[2], t) == -μ_nd * r[2] / (r[1]^2 + r[2]^2 + r[3]^2)^(3 / 2)
                          + (T_max_nd / m) * δ * sin(θ) * sin(φ)
        )
    else
        @constraint(model,
            ∂(v[3], t) == -μ_nd * r[3] / (r[1]^2 + r[2]^2 + r[3]^2)^(3 / 2)
                          + (T_max_nd / m) * δ * cos(θ)
        )
    end
end

# Dynamics: dm/dt
@constraint(model, ∂(m, t) == -(T_max_nd / c_nd) * δ)

# Path constraint: mass must stay positive
@constraint(model, [t = t], m >= m_min_nd)

println("Model built successfully!")
println("  State variables: 7 (r[3], v[3], m)")
println("  Control variables: 3 (θ, φ, δ)")
println("  Collocation nodes per element: 4")

# =============================================================================
# SOLVE
# =============================================================================
println("\n" * "="^70)
println("Solving...")
println("="^70)

start_time = time()
optimize!(model)
solve_time = time() - start_time

# =============================================================================
# EXTRACT RESULTS
# =============================================================================
success = termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
obj_value = objective_value(model)

println("\n" * "="^70)
println("Optimization Results")
println("="^70)
println("  Success: $success")
println("  Termination status: $(termination_status(model))")
println("  Solve time: $(round(solve_time, digits=2)) s")
println("  Objective value: $(round(obj_value, digits=6))")

# Get solution
t_vals = value(t)
r_vals = value.(r)
v_vals = value.(v)
m_vals = value(m)
θ_vals = value(θ)
φ_vals = value(φ)
δ_vals = value(δ)

N_actual = length(t_vals)
X_sol = zeros(7, N_actual)
U_sol = zeros(4, N_actual)

for i in 1:N_actual
    X_sol[1:3, i] = [r_vals[j][i] for j in 1:3]
    X_sol[4:6, i] = [v_vals[j][i] for j in 1:3]
    X_sol[7, i] = m_vals[i]

    U_sol[1, i] = θ_vals[i]
    U_sol[2, i] = φ_vals[i]
    U_sol[3, i] = δ_vals[i]
    U_sol[4, i] = T_max_nd * δ_vals[i]
end

# Re-dimensionalize
X_sol[1:3, :] .*= L_scale
X_sol[4:6, :] .*= V_scale
X_sol[7, :] .*= M_scale
t_vals_dim = t_vals .* T_scale

# Print key results
println("\n" * "="^70)
println("Solution Summary")
println("="^70)
println("  Initial mass: $(round(X_sol[7, 1], digits=2)) kg")
println("  Final mass: $(round(X_sol[7, end], digits=2)) kg")
println("  Mass used: $(round(X_sol[7, 1] - X_sol[7, end], digits=2)) kg")
println("")
println("  Initial position: [$(round(X_sol[1, 1]/AU_m, digits=4)), $(round(X_sol[2, 1]/AU_m, digits=4)), $(round(X_sol[3, 1]/AU_m, digits=4))] AU")
println("  Final position: [$(round(X_sol[1, end]/AU_m, digits=4)), $(round(X_sol[2, end]/AU_m, digits=4)), $(round(X_sol[3, end]/AU_m, digits=4))] AU")
println("")
println("  Initial velocity: [$(round(X_sol[4, 1]/1000, digits=2)), $(round(X_sol[5, 1]/1000, digits=2)), $(round(X_sol[6, 1]/1000, digits=2))] km/s")
println("  Final velocity: [$(round(X_sol[4, end]/1000, digits=2)), $(round(X_sol[5, end]/1000, digits=2)), $(round(X_sol[6, end]/1000, digits=2))] km/s")

# Check constraint violations
println("\n" * "="^70)
println("Constraint Violations")
println("="^70)

# Position error
pos_error = norm(X_sol[1:3, end] - xf[1:3])
println("  Final position error: $(round(pos_error / 1000, digits=3)) km")

# Mass constraint
min_mass = minimum(X_sol[7, :])
println("  Minimum mass: $(round(min_mass, digits=2)) kg")
if min_mass < 500.0
    println("  ⚠️  Mass constraint violated!")
else
    println("  ✓ Mass constraint satisfied")
end

# Throttle stats
avg_throttle = mean(U_sol[3, :])
max_throttle = maximum(U_sol[3, :])
println("\n" * "="^70)
println("Control Statistics")
println("="^70)
println("  Average throttle: $(round(avg_throttle, digits=3))")
println("  Max throttle: $(round(max_throttle, digits=3))")

# Thrust direction stats
thrust_mags = [norm([sin(U_sol[1, i]) * cos(U_sol[2, i]),
                     sin(U_sol[1, i]) * sin(U_sol[2, i]),
                     cos(U_sol[1, i])]) for i in 1:N_actual]
avg_thrust_mag = mean(thrust_mags)
max_thrust_deviation = maximum(abs.(thrust_mags .- 1.0))
println("  Thrust direction normalization error: $(round(max_thrust_deviation, digits=10))")
if max_thrust_deviation < 1e-6
    println("  ✓ Thrust direction properly normalized")
end

println("\n" * "="^70)
println("Test Complete!")
println("="^70)

if success && pos_error < 1e6 && min_mass >= 500.0
    println("✅ ALL CHECKS PASSED")
    println("InfiniteOpt low-thrust optimizer is working correctly!")
else
    println("⚠️  Some checks failed - review results above")
end

println("\nNext steps:")
println("  1. Run full GTOC12 example with hit_mine_an_asteroid.jl")
println("  2. Compare with shooting method results")
println("  3. Implement MEE coordinate version")
println("  4. Validate with GTOC12 verification program")
println("="^70)
