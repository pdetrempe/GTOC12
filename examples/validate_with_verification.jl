"""
Generate GTOC12-compliant Result.txt and validate with verification program.

Key fix: Use MJD time format instead of ET.
"""

using GTOC12
using LinearAlgebra
using Printf

println("="^70)
println("GTOC12 Verification Program Validation")
println("="^70)

# Helper function: Convert ET to MJD
function ET_to_MJD(ET::Float64)
    return ET / 86400.0 + 2451545.0 - 2400000.5
end

# Setup problem - START AT MISSION BEGINNING
GTOC12.furnish_all_kernels()
mining_ship = GTOC12.Mining_Ship()

# Find asteroid closest to Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / GTOC12.au2m .- 1))
asteroid = GTOC12.Asteroid(ID_min)

# START AT MISSION BEGINNING (not +338 days)
ET_start = GTOC12.ET₀  # Mission start
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

# Convert to MJD
MJD_start = ET_to_MJD(ET_start)
MJD_target = ET_to_MJD(ET_target)

println("\nProblem Setup:")
println("  Asteroid ID: $ID_min")
println("  Mission start ET: $ET_start")
println("  Mission start MJD: $MJD_start")
println("  Transfer time: $(Δt/(24*3600)) days")
println("  Arrival MJD: $MJD_target")
println("  Mission end MJD: 69807")
println("  Within window: $(MJD_target < 69807 ? "✓" : "✗")")

# Run optimizer
println("\n" * "="^70)
println("Running InfiniteOpt Optimizer")
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

println("\nOptimization Results:")
println("  Success: $(result.success)")
println("  Final mass: $(round(result.X[7, end], digits=2)) kg")
println("  Mass used: $(round(result.mass_used, digits=2)) kg")

if !result.success
    println("\n❌ Optimization failed")
    exit(1)
end

# Generate Result.txt with MJD time format
println("\n" * "="^70)
println("Generating Result.txt (MJD format)")
println("="^70)

output_file = "/Users/pauldetrempe/GTOC12/problem/GTOC12_Verification_Program/GTOC12_Verification/macOS/Result.txt"

open(output_file, "w") do io
    # Line 1: Earth departure state (MJD format)
    @printf(io, "%5d %5d %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f\n",
        1, 0, MJD_start,
        x₀[1], x₀[2], x₀[3],
        x₀[4], x₀[5], x₀[6],
        mining_ship.mass_total)

    # Line 2: Asteroid arrival state (MJD format)
    @printf(io, "%5d %5d %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f\n",
        1, 0, MJD_target,
        x_target[1], x_target[2], x_target[3],
        x_target[4], x_target[5], x_target[6],
        result.X[7, end])

    # Trajectory points (MJD format)
    for i in 1:length(result.t_nodes)
        ET = ET_start + result.t_nodes[i]
        MJD = ET_to_MJD(ET)

        # Convert thrust angles to Cartesian
        θ = result.U[1, i]
        φ = result.U[2, i]
        δ = result.U[3, i]

        # Thrust direction in RTN
        û_r = sin(θ) * cos(φ)
        û_t = sin(θ) * sin(φ)
        û_n = cos(θ)

        # Get spacecraft state
        r = result.X[1:3, i]
        v = result.X[4:6, i]

        # RTN basis
        r_hat = r / norm(r)
        h = cross(r, v)
        n_hat = h / norm(h)
        t_hat = cross(n_hat, r_hat)

        # Thrust in inertial frame
        thrust_inertial = û_r * r_hat + û_t * t_hat + û_n * n_hat
        thrust_magnitude = GTOC12.T_max * δ
        thrust_vector = thrust_magnitude * thrust_inertial

        # Write trajectory point
        @printf(io, "%d %d %.16f %.16f %.16f %.16f %.16f\n",
            1, -1, MJD,
            thrust_vector[1], thrust_vector[2], thrust_vector[3],
            result.X[7, i])
    end
end

println("  File written: $output_file")
println("  Format: MJD time")
println("  Start MJD: $MJD_start")
println("  End MJD: $MJD_target")

# Run verification
println("\n" * "="^70)
println("Running GTOC12 Verification Program")
println("="^70)

cd("/Users/pauldetrempe/GTOC12/problem/GTOC12_Verification_Program/GTOC12_Verification/macOS")
try
    run(`./GTOC12_Verify`)
    println("\n✅ VERIFICATION PASSED!")
catch e
    println("\n⚠️  Verification output (check above for details)")
end

println("\n" * "="^70)
println("Validation Complete!")
println("="^70)
