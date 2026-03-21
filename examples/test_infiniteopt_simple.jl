"""
Simple test of InfiniteOpt low-thrust optimizer

This script tests the InfiniteOpt implementation with a simple intercept problem:
- Start from Earth
- Intercept an asteroid (position only, velocity free)
- Compare with shooting method from hit_mine_an_asteroid.jl
"""

using GTOC12
using LinearAlgebra
using Plots

println("="^70)
println("Testing InfiniteOpt Low-Thrust Optimizer")
println("="^70)

# Load InfiniteOpt module
using GTOC12.InfiniteOptLowThrust

# Initialize Mining Ship
mining_ship = GTOC12.Mining_Ship()
println("Mining ship initial mass: $(mining_ship.mass_total) kg")

# Find asteroid closest (in terms of orbital energy) to the Earth
furnish_all_kernels()
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1))
asteroid = Asteroid(ID_min)
println("Target asteroid ID: $ID_min")

# Initial conditions (from Earth)
ET_start = GTOC12.ET₀ + 338*24*3600
x_Earth = get_body_state(GTOC12.Earth; ET=ET_start)
println("Earth position: $(x_Earth[1:3] / 1e9) km (×10⁹)")

# Launch with some ΔV (from hit_mine_an_asteroid.jl)
DV₀ = [159.31762519335462; 1505.6951917938995; 3569.929588119754]
x₀ = x_Earth + [0; 0; 0; DV₀[:]]
println("Launch ΔV: $(norm(DV₀)) m/s")

# Transfer time
Δt = 0.35 * 365*24*3600  # seconds
println("Transfer time: $(Δt / (365*24*3600)) years")

# Target state
ET_target = ET_start + Δt
x_target = get_body_state(asteroid; ET=ET_target)
println("Target position: $(x_target[1:3] / 1e9) km (×10⁹)")

# Prepare initial state for InfiniteOpt (with mass)
x0_with_mass = vcat(x₀, mining_ship.mass_total)

# Run InfiniteOpt optimizer (intercept only)
println("\n" * "="^70)
println("Running InfiniteOpt Optimizer (Intercept Only)")
println("="^70)

N = 50  # discretization points
result = infiniteopt_low_thrust_cartesian(
    x0_with_mass,
    x_target,
    Δt,
    N,
    GTOC12.μ_☉;
    T_max=GTOC12.T_max,
    Isp=Float64(GTOC12.Isp),
    g0=GTOC12.g0,
    m_min=500.0,
    intercept_only=true,
    max_iter=3000,
    verbose=true
)

println("\n" * "="^70)
println("Optimization Results:")
println("="^70)
println("  Success: $(result.success)")
println("  Solve time: $(round(result.solve_time, digits=2)) s")
println("  Objective value: $(round(result.objective_value, digits=6))")
println("  Initial mass: $(round(result.X[7, 1], digits=2)) kg")
println("  Final mass: $(round(result.X[7, end], digits=2)) kg")
println("  Mass used: $(round(result.mass_used, digits=2)) kg")
println("")
println("  Final position error: $(round(result.constraint_breakdown.position / 1e3, digits=3)) km")
println("  Final velocity error: $(round(result.constraint_breakdown.velocity, digits=3)) m/s")
println("  Max constraint violation: $(round(result.constraint_violation / 1e3, digits=6)) km")
println("="^70)

# Extract trajectory
r_traj = result.X[1:3, :] / 1e9  # km (×10⁹)
v_traj = result.X[4:6, :]
m_traj = result.X[7, :]
t_traj = result.t_nodes / (365*24*3600)  # years

# Extract controls (convert from spherical to Cartesian)
θ_vals = result.U[1, :]
φ_vals = result.U[2, :]
δ_vals = result.U[3, :]

thrust_directions = zeros(3, length(θ_vals))
for i in 1:length(θ_vals)
    thrust_directions[1, i] = sin(θ_vals[i]) * cos(φ_vals[i])
    thrust_directions[2, i] = sin(θ_vals[i]) * sin(φ_vals[i])
    thrust_directions[3, i] = cos(θ_vals[i])
end

# Plotting
println("\nGenerating plots...")

# Plot 1: 3D trajectory
p1 = plot(title="3D Trajectory", xlabel="X (km ×10⁹)", ylabel="Y (km ×10⁹)", zlabel="Z (km ×10⁹)",
    legend=:topright, aspect_ratio=:equal)

# Earth orbit
ET_plot = range(ET_start, ET_target, length=100)
earth_traj = hcat([get_body_state(GTOC12.Earth; ET=et)[1:3] for et in ET_plot]...) / 1e9
plot!(p1, earth_traj[1, :], earth_traj[2, :], earth_traj[3, :], label="Earth orbit", color=:blue, linewidth=1)

# Asteroid orbit
asteroid_traj = hcat([get_body_state(asteroid; ET=et)[1:3] for et in ET_plot]...) / 1e9
plot!(p1, asteroid_traj[1, :], asteroid_traj[2, :], asteroid_traj[3, :], label="Asteroid orbit", color=:red, linewidth=1)

# Spacecraft trajectory
plot!(p1, r_traj[1, :], r_traj[2, :], r_traj[3, :], label="Spacecraft", color=:green, linewidth=2)
scatter!(p1, [r_traj[1, 1]], [r_traj[2, 1]], [r_traj[3, 1]], label="Start", color=:blue, markersize=5)
scatter!(p1, [r_traj[1, end]], [r_traj[2, end]], [r_traj[3, end]], label="End", color=:red, markersize=5)

# Plot 2: Mass history
p2 = plot(title="Mass History", xlabel="Time (years)", ylabel="Mass (kg)", legend=:topright)
plot!(p2, t_traj, m_traj, label="Spacecraft mass", color=:blue, linewidth=2)

# Plot 3: Throttle history
p3 = plot(title="Throttle History", xlabel="Time (years)", ylabel="Throttle δ", legend=:topright, ylim=(0, 1.1))
plot!(p3, t_traj, δ_vals, label="Throttle", color=:orange, linewidth=2)

# Plot 4: Thrust magnitude
p4 = plot(title="Thrust Magnitude", xlabel="Time (years)", ylabel="Thrust (N)", legend=:topright)
thrust_mag = GTOC12.T_max .* δ_vals
plot!(p4, t_traj, thrust_mag, label="Thrust", color=:red, linewidth=2)

# Combine plots
plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 900))
savefig("test_infiniteopt_simple_results.png")
println("Plots saved to test_infiniteopt_simple_results.png")

println("\n" * "="^70)
println("Test complete!")
println("="^70)
