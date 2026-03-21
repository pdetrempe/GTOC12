"""
Fix the Result.txt to use MJD format instead of ET format.

This takes the existing working solution and just converts the time format.
"""

using GTOC12
using LinearAlgebra
using Printf

# Helper function: Convert ET to MJD
function ET_to_MJD(ET::Float64)
    return ET / 86400.0 + 2451545.0 - 2400000.5
end

# Use the working test problem parameters
GTOC12.furnish_all_kernels()
mining_ship = GTOC12.Mining_Ship()

# Find asteroid
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / GTOC12.au2m .- 1))
asteroid = GTOC12.Asteroid(ID_min)

# Use the working start time
ET_start = GTOC12.ET₀ + 338*24*3600  # This worked in our tests
x_Earth = GTOC12.get_body_state(GTOC12.Earth; ET=ET_start)

DV₀ = [159.31762519335462; 1505.6951917938995; 3569.929588119754]
x₀ = x_Earth + [0; 0; 0; DV₀[:]]

Δt = 0.35 * 365*24*3600
ET_target = ET_start + Δt
x_target = GTOC12.get_body_state(asteroid; ET=ET_target)
x0_with_mass = vcat(x₀, mining_ship.mass_total)

# Convert to MJD
MJD_start = ET_to_MJD(ET_start)
MJD_target = ET_to_MJD(ET_target)

println("="^70)
println("Fixing Result.txt Format: ET → MJD")
println("="^70)
println("\nTime Conversion:")
println("  ET start:  $ET_start")
println("  MJD start: $MJD_start")
println("  MJD end:   $MJD_target")
println("  Mission window: 64328 - 69807 MJD")

# Check if within mission window
if MJD_start < 64328 || MJD_target > 69807
    println("\n⚠️  WARNING: Times outside mission window!")
    println("  Adjusting start time to mission beginning...")

    # Use mission start instead
    ET_start = GTOC12.ET₀
    MJD_start = ET_to_MJD(ET_start)
    MJD_target = ET_to_MJD(ET_start + Δt)

    # Recalculate states
    x_Earth = GTOC12.get_body_state(GTOC12.Earth; ET=ET_start)
    x₀ = x_Earth + [0; 0; 0; DV₀[:]]
    x0_with_mass = vcat(x₀, mining_ship.mass_total)
    ET_target = ET_start + Δt
    x_target = GTOC12.get_body_state(asteroid; ET=ET_target)

    println("  New MJD start: $MJD_start")
    println("  New MJD end:   $MJD_target")
end

println("\n✓ Times are within mission window")

# Run optimizer
println("\n" * "="^70)
println("Running Optimizer")
println("="^70)

N = 50
result = GTOC12.InfiniteOptLowThrust.infiniteopt_low_thrust_cartesian(
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
    verbose=false
)

println("  Success: $(result.success)")
println("  Final mass: $(round(result.X[7, end], digits=2)) kg")

if !result.success
    println("\n❌ Failed - cannot generate Result.txt")
    exit(1)
end

# Generate Result.txt
println("\n" * "="^70)
println("Generating Result.txt")
println("="^70)

output_file = "/Users/pauldetrempe/GTOC12/problem/GTOC12_Verification_Program/GTOC12_Verification/macOS/Result.txt"

open(output_file, "w") do io
    # Header lines with MJD
    # Line 1: Earth departure (asteroid ID = 0 for Earth)
    @printf(io, "%5d %5d %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f\n",
        1, 0, MJD_start,
        x₀[1], x₀[2], x₀[3],
        x₀[4], x₀[5], x₀[6],
        mining_ship.mass_total)

    # Line 2: Asteroid arrival (use actual asteroid ID)
    @printf(io, "%5d %5d %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f\n",
        1, ID_min, MJD_target,
        x_target[1], x_target[2], x_target[3],
        x_target[4], x_target[5], x_target[6],
        result.X[7, end])

    # Trajectory points with MJD (skip first point if it's at t=0, already in header)
    start_idx = (result.t_nodes[1] ≈ 0.0) ? 2 : 1
    for i in start_idx:length(result.t_nodes)
        MJD = ET_to_MJD(ET_start + result.t_nodes[i])

        θ = result.U[1, i]
        φ = result.U[2, i]
        δ = result.U[3, i]

        û_r = sin(θ) * cos(φ)
        û_t = sin(θ) * sin(φ)
        û_n = cos(θ)

        r = result.X[1:3, i]
        v = result.X[4:6, i]

        r_hat = r / norm(r)
        h = cross(r, v)
        n_hat = h / norm(h)
        t_hat = cross(n_hat, r_hat)

        thrust_inertial = û_r * r_hat + û_t * t_hat + û_n * n_hat
        thrust_vector = GTOC12.T_max * δ * thrust_inertial

        # Trajectory point during cruise (asteroid_id = -1 for control lines)
        # CRITICAL: Control lines have 6 fields only (NO mass field!)
        @printf(io, "%5d %5d %25.12f %25.16f %25.16f %25.16f\n",
            1, -1, MJD,
            thrust_vector[1], thrust_vector[2], thrust_vector[3])
    end
end

println("  ✓ File written: $output_file")
println("  Format: MJD (not ET)")

# Verify
println("\n" * "="^70)
println("Running GTOC12 Verification")
println("="^70)

cd("/Users/pauldetrempe/GTOC12/problem/GTOC12_Verification_Program/GTOC12_Verification/macOS")
run(`./GTOC12_Verify`)

println("\n" * "="^70)
println("Complete!")
println("="^70)
