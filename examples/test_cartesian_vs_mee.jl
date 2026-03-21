"""
Test Cartesian vs MEE Coordinate Formulations

Compares InfiniteOpt low-thrust optimization using:
1. Cartesian coordinates (r, v, m)
2. Modified Equinoctial Elements (p, f, g, h, k, L, m)

Evaluates numerical conditioning and convergence characteristics.
"""

using GTOC12
using Plots
using Printf
using LinearAlgebra

println("="^70)
println("Cartesian vs MEE Comparison Test")
println("="^70)

# Use GTOC12 infrastructure
GTOC12.furnish_all_kernels()

# Setup problem (same as test_infiniteopt_simple.jl)
mining_ship = GTOC12.Mining_Ship()
println("\nMining ship initial mass: $(mining_ship.mass_total) kg")

# Find asteroid closest to Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / GTOC12.au2m .- 1))
asteroid = GTOC12.Asteroid(ID_min)
println("Target asteroid ID: $ID_min")

# Initial conditions (from Earth)
ET_start = GTOC12.ET₀ + 338*24*3600
x_Earth = GTOC12.get_body_state(GTOC12.Earth; ET=ET_start)

# Launch with some ΔV
DV₀ = [159.31762519335462; 1505.6951917938995; 3569.929588119754]
x₀ = x_Earth + [0; 0; 0; DV₀[:]]

# Transfer time
Δt = 0.35 * 365*24*3600  # seconds

# Target state
ET_target = ET_start + Δt
x_target = GTOC12.get_body_state(asteroid; ET=ET_target)

# Prepare states
x0_with_mass = vcat(x₀, mining_ship.mass_total)

println("\nEarth position: $([round(x_Earth[i]/1e9, digits=2) for i in 1:3]) km (×10⁹)")
println("Launch ΔV: $(norm(DV₀)) m/s")
println("Transfer time: $(Δt / (365*24*3600)) years")
println("Target position: $([round(x_target[i]/1e9, digits=2) for i in 1:3]) km (×10⁹)")

# Common parameters
N = 50
params = (
    T_max=GTOC12.T_max,
    Isp=Float64(GTOC12.Isp),
    g0=GTOC12.g0,
    m_min=500.0,
    intercept_only=true,
    max_iter=3000,
    verbose=false
)

println("\n" * "="^70)
println("Test 1: Cartesian Coordinate Formulation")
println("="^70)

result_cart = GTOC12.InfiniteOptLowThrust.infiniteopt_low_thrust_cartesian(
    x0_with_mass,
    x_target,
    Δt,
    N,
    GTOC12.μ_☉;
    params...
)

println("\nCartesian Results:")
println("  Success: $(result_cart.success)")
println("  Solve time: $(round(result_cart.solve_time, digits=2)) s")
println("  Objective value: $(round(result_cart.objective_value, digits=6))")
println("  Final mass: $(round(result_cart.X[7, end], digits=2)) kg")
println("  Mass used: $(round(result_cart.mass_used, digits=2)) kg")
println("  Constraint violation: $(round(result_cart.constraint_violation, digits=10))")
println("    - Position: $(round(result_cart.constraint_breakdown.position, digits=3)) m")
println("    - Velocity: $(round(result_cart.constraint_breakdown.velocity, digits=3)) m/s")

println("\n" * "="^70)
println("Test 2: MEE Coordinate Formulation")
println("="^70)

result_mee = GTOC12.InfiniteOptLowThrust.infiniteopt_low_thrust_mee(
    x0_with_mass,
    x_target,
    Δt,
    N,
    GTOC12.μ_☉;
    params...,
    convert_to_mee=true
)

println("\nMEE Results:")
println("  Success: $(result_mee.success)")
println("  Solve time: $(round(result_mee.solve_time, digits=2)) s")
println("  Objective value: $(round(result_mee.objective_value, digits=6))")
println("  Final mass: $(round(result_mee.X[7, end], digits=2)) kg")
println("  Mass used: $(round(result_mee.mass_used, digits=2)) kg")
println("  Constraint violation: $(round(result_mee.constraint_violation, digits=10))")
println("    - MEE elements: $(round(result_mee.constraint_breakdown.position, digits=10))")
println("    - True longitude: $(round(result_mee.constraint_breakdown.velocity, digits=10))")

println("\n" * "="^70)
println("Comparison Analysis")
println("="^70)

# Convert MEE to Cartesian for comparison
X_mee_cart = zeros(7, size(result_mee.X, 2))
for i in 1:size(result_mee.X, 2)
    mee_state = result_mee.X[1:6, i]
    cart_state = GTOC12.MEE2Cartesian(mee_state; μ=GTOC12.μ_☉)
    X_mee_cart[1:6, i] = cart_state
    X_mee_cart[7, i] = result_mee.X[7, i]
end

# Compare final states
println("\nFinal State Comparison:")
println("  Position difference: $(norm(result_cart.X[1:3, end] - X_mee_cart[1:3, end])/1000) km")
println("  Velocity difference: $(norm(result_cart.X[4:6, end] - X_mee_cart[4:6, end])) m/s")
println("  Mass difference: $(abs(result_cart.X[7, end] - X_mee_cart[7, end])) kg")

# Compare mass usage
mass_diff = abs(result_cart.mass_used - result_mee.mass_used)
mass_diff_pct = 100 * mass_diff / result_cart.mass_used
println("\nMass Usage:")
println("  Cartesian: $(round(result_cart.mass_used, digits=2)) kg")
println("  MEE:       $(round(result_mee.mass_used, digits=2)) kg")
println("  Difference: $(round(mass_diff, digits=2)) kg ($(round(mass_diff_pct, digits=2))%)")

# Compare solve times
time_ratio = result_mee.solve_time / result_cart.solve_time
println("\nSolve Time:")
println("  Cartesian: $(round(result_cart.solve_time, digits=2)) s")
println("  MEE:       $(round(result_mee.solve_time, digits=2)) s")
println("  Ratio (MEE/Cart): $(round(time_ratio, digits=2))x")

# Numerical conditioning assessment
println("\nNumerical Conditioning:")
println("  Cartesian constraint violation: $(round(result_cart.constraint_violation, digits=12))")
println("  MEE constraint violation:       $(round(result_mee.constraint_violation, digits=12))")

if result_cart.constraint_violation < result_mee.constraint_violation
    better = "Cartesian"
    ratio = result_mee.constraint_violation / result_cart.constraint_violation
else
    better = "MEE"
    ratio = result_cart.constraint_violation / result_mee.constraint_violation
end
println("  Better conditioning: $better ($(round(ratio, digits=2))x better)")

# Average throttle comparison
avg_throttle_cart = sum(result_cart.U[3, :]) / size(result_cart.U, 2)
avg_throttle_mee = sum(result_mee.U[3, :]) / size(result_mee.U, 2)
println("\nControl Statistics:")
println("  Cartesian avg throttle: $(round(avg_throttle_cart, digits=3))")
println("  MEE avg throttle:       $(round(avg_throttle_mee, digits=3))")

println("\n" * "="^70)
println("Plotting Results")
println("="^70)

# Plot trajectories
p1 = plot(result_cart.X[1, :]/1e9, result_cart.X[2, :]/1e9,
    label="Cartesian", linewidth=2, xlabel="X (km ×10⁹)", ylabel="Y (km ×10⁹)",
    title="Trajectory Comparison (XY Plane)")
plot!(p1, X_mee_cart[1, :]/1e9, X_mee_cart[2, :]/1e9,
    label="MEE", linewidth=2, linestyle=:dash)
scatter!(p1, [x_Earth[1]/1e9], [x_Earth[2]/1e9], label="Earth", markersize=6, color=:blue)
scatter!(p1, [x_target[1]/1e9], [x_target[2]/1e9], label="Target", markersize=6, color=:red)

# Plot mass histories
p2 = plot(result_cart.t_nodes/(24*3600), result_cart.X[7, :],
    label="Cartesian", linewidth=2, xlabel="Time (days)", ylabel="Mass (kg)",
    title="Mass History")
plot!(p2, result_mee.t_nodes/(24*3600), result_mee.X[7, :],
    label="MEE", linewidth=2, linestyle=:dash)

# Plot throttle histories
p3 = plot(result_cart.t_nodes/(24*3600), result_cart.U[3, :],
    label="Cartesian", linewidth=2, xlabel="Time (days)", ylabel="Throttle",
    title="Throttle History")
plot!(p3, result_mee.t_nodes/(24*3600), result_mee.U[3, :],
    label="MEE", linewidth=2, linestyle=:dash)

# Plot position error over time
pos_error = zeros(size(result_cart.X, 2))
for i in 1:size(result_cart.X, 2)
    pos_error[i] = norm(result_cart.X[1:3, i] - X_mee_cart[1:3, i])
end
p4 = plot(result_cart.t_nodes/(24*3600), pos_error/1000,
    label="Position Error", linewidth=2, xlabel="Time (days)", ylabel="Error (km)",
    title="Cartesian vs MEE Position Difference")

# Combine plots
plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 900))
savefig("cartesian_vs_mee_comparison.png")
println("Plots saved to cartesian_vs_mee_comparison.png")

println("\n" * "="^70)
println("Summary and Recommendations")
println("="^70)

if result_cart.success && result_mee.success
    println("✅ Both formulations converged successfully")

    if mass_diff_pct < 1.0
        println("✅ Mass usage agrees within $(round(mass_diff_pct, digits=2))% - solutions are consistent")
    else
        println("⚠️  Mass usage differs by $(round(mass_diff_pct, digits=2))% - may indicate local minima")
    end

    if better == "MEE"
        println("✅ MEE has better numerical conditioning ($(round(ratio, digits=2))x)")
        println("   Recommendation: Use MEE for production runs")
    else
        println("ℹ️  Cartesian has better numerical conditioning ($(round(ratio, digits=2))x)")
        println("   Recommendation: Both formulations are viable")
    end

    if time_ratio < 1.2
        println("✅ Solve times are comparable (within 20%)")
    elseif time_ratio > 1.5
        println("⚠️  MEE is slower by $(round((time_ratio-1)*100, digits=0))%")
    else
        println("ℹ️  MEE is slightly slower by $(round((time_ratio-1)*100, digits=0))%")
    end

elseif result_cart.success && !result_mee.success
    println("⚠️  Cartesian converged but MEE did not")
    println("   Recommendation: Use Cartesian formulation")
elseif !result_cart.success && result_mee.success
    println("⚠️  MEE converged but Cartesian did not")
    println("   Recommendation: Use MEE formulation")
else
    println("❌ Both formulations failed to converge")
    println("   Recommendation: Adjust problem parameters or initial guess")
end

println("\n" * "="^70)
println("Test Complete!")
println("="^70)
