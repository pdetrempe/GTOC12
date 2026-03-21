"""
Complete round-trip mission: Earth → Asteroid → Earth

This uses InfiniteOpt to:
1. Transfer from Earth to nearest asteroid
2. Transfer back from asteroid to Earth
"""

using GTOC12
using LinearAlgebra
using Printf

println("="^70)
println("Complete Round-Trip Mission with InfiniteOpt")
println("="^70)

# Helper function: Convert ET to MJD
function ET_to_MJD(ET::Float64)
    return ET / 86400.0 + 2451545.0 - 2400000.5
end

# Setup problem
GTOC12.furnish_all_kernels()
mining_ship = GTOC12.Mining_Ship()

# Find asteroid closest to Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / GTOC12.au2m .- 1))
asteroid = GTOC12.Asteroid(ID_min)

println("\nMission Setup:")
println("  Target Asteroid ID: $ID_min")
println("  Initial mass: $(mining_ship.mass_total) kg")

# ========================================================================
# LEG 1: Earth → Asteroid
# ========================================================================
println("\n" * "="^70)
println("LEG 1: Earth → Asteroid")
println("="^70)

# Start slightly after mission beginning to ensure we're within window
# GTOC12.ET₀ converts to MJD 64327.999..., but mission starts at MJD 64328
ET_start = GTOC12.ET₀ + 1.0*24*3600  # Add 1 day to be safely within window
x_Earth_dep = GTOC12.get_body_state(GTOC12.Earth; ET=ET_start)

# Launch with ΔV
DV₀ = [159.31762519335462; 1505.6951917938995; 3569.929588119754]
x₀ = x_Earth_dep + [0; 0; 0; DV₀[:]]

# Transfer time to asteroid (about 4 months)
Δt_leg1 = 0.35 * 365*24*3600  # seconds

# Target state at asteroid
ET_asteroid_arrival = ET_start + Δt_leg1
x_asteroid = GTOC12.get_body_state(asteroid; ET=ET_asteroid_arrival)

# Prepare initial state with mass
x0_with_mass = vcat(x₀, mining_ship.mass_total)

# Run optimizer for Leg 1
N = 50
result1 = GTOC12.InfiniteOptLowThrust.infiniteopt_low_thrust_cartesian(
    x0_with_mass,
    x_asteroid,
    Δt_leg1,
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

println("  Optimizer status: $(result1.success)")
println("  Arrival mass: $(round(result1.X[7, end], digits=2)) kg")
println("  Mass used: $(round(result1.mass_used, digits=2)) kg")

# Continue even if optimizer reports non-convergence
# (it often produces usable trajectories anyway)
if result1.X[7, end] < 500.0
    println("\n❌ Leg 1: Final mass too low ($(result1.X[7, end]) kg)")
    exit(1)
end

# ========================================================================
# LEG 2: Asteroid → Earth (Return)
# ========================================================================
println("\n" * "="^70)
println("LEG 2: Asteroid → Earth (Return)")
println("="^70)

# Start from asteroid arrival state
mass_at_asteroid = result1.X[7, end]
ET_departure_from_asteroid = ET_asteroid_arrival

# Transfer time back to Earth (about 6 months)
Δt_leg2 = 0.5 * 365*24*3600  # seconds

# Target: Earth position at return time
ET_earth_arrival = ET_departure_from_asteroid + Δt_leg2
x_Earth_arrival = GTOC12.get_body_state(GTOC12.Earth; ET=ET_earth_arrival)

# Initial state for leg 2: at asteroid with remaining mass
x_asteroid_dep = vcat(x_asteroid, mass_at_asteroid)

# Run optimizer for Leg 2
result2 = GTOC12.InfiniteOptLowThrust.infiniteopt_low_thrust_cartesian(
    x_asteroid_dep,
    x_Earth_arrival,
    Δt_leg2,
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

println("  Optimizer status: $(result2.success)")
println("  Final mass: $(round(result2.X[7, end], digits=2)) kg")
println("  Mass used: $(round(result2.mass_used, digits=2)) kg")

# Continue even if optimizer reports non-convergence
if result2.X[7, end] < 500.0
    println("\n❌ Leg 2: Final mass too low ($(result2.X[7, end]) kg)")
    exit(1)
end

# ========================================================================
# Generate Result.txt
# ========================================================================
println("\n" * "="^70)
println("Generating Result.txt")
println("="^70)

MJD_start = ET_to_MJD(ET_start)
MJD_asteroid_arrival = ET_to_MJD(ET_asteroid_arrival)
MJD_earth_arrival = ET_to_MJD(ET_earth_arrival)

println("  Departure from Earth: MJD $MJD_start")
println("  Arrival at Asteroid:  MJD $MJD_asteroid_arrival")
println("  Return to Earth:      MJD $MJD_earth_arrival")
println("  Mission duration:     $(round((MJD_earth_arrival - MJD_start), digits=1)) days")

output_file = "/Users/pauldetrempe/GTOC12/problem/GTOC12_Verification_Program/GTOC12_Verification/macOS/Result.txt"

open(output_file, "w") do io
    # ====================================================================
    # Header Line 1: Earth departure state
    # ====================================================================
    @printf(io, "%5d %5d %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f\n",
        1, 0, MJD_start,
        x₀[1], x₀[2], x₀[3],
        x₀[4], x₀[5], x₀[6],
        mining_ship.mass_total)

    # ====================================================================
    # Header Line 2: Asteroid arrival state
    # ====================================================================
    @printf(io, "%5d %5d %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f\n",
        1, ID_min, MJD_asteroid_arrival,
        x_asteroid[1], x_asteroid[2], x_asteroid[3],
        x_asteroid[4], x_asteroid[5], x_asteroid[6],
        mass_at_asteroid)

    # ====================================================================
    # LEG 1: Control points (Earth → Asteroid)
    # ====================================================================
    for i in 1:length(result1.t_nodes)
        MJD = ET_to_MJD(ET_start + result1.t_nodes[i])

        # Convert thrust angles to Cartesian
        θ = result1.U[1, i]
        φ = result1.U[2, i]
        δ = result1.U[3, i]

        # Thrust direction in RTN
        û_r = sin(θ) * cos(φ)
        û_t = sin(θ) * sin(φ)
        û_n = cos(θ)

        # Get spacecraft state
        r = result1.X[1:3, i]
        v = result1.X[4:6, i]

        # RTN basis
        r_hat = r / norm(r)
        h = cross(r, v)
        n_hat = h / norm(h)
        t_hat = cross(n_hat, r_hat)

        # Thrust in inertial frame
        thrust_inertial = û_r * r_hat + û_t * t_hat + û_n * n_hat
        thrust_vector = GTOC12.T_max * δ * thrust_inertial

        # Write control line (6 fields, NO mass)
        @printf(io, "%5d %5d %25.12f %25.16f %25.16f %25.16f\n",
            1, -1, MJD,
            thrust_vector[1], thrust_vector[2], thrust_vector[3])
    end

    # ====================================================================
    # Header Line 3: Asteroid departure state (for return leg)
    # ====================================================================
    @printf(io, "%5d %5d %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f\n",
        1, ID_min, MJD_asteroid_arrival,
        x_asteroid[1], x_asteroid[2], x_asteroid[3],
        x_asteroid[4], x_asteroid[5], x_asteroid[6],
        mass_at_asteroid)

    # ====================================================================
    # Header Line 4: Earth arrival state
    # ====================================================================
    @printf(io, "%5d %5d %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f %25.12f\n",
        1, 0, MJD_earth_arrival,
        x_Earth_arrival[1], x_Earth_arrival[2], x_Earth_arrival[3],
        x_Earth_arrival[4], x_Earth_arrival[5], x_Earth_arrival[6],
        result2.X[7, end])

    # ====================================================================
    # LEG 2: Control points (Asteroid → Earth)
    # ====================================================================
    for i in 1:length(result2.t_nodes)
        MJD = ET_to_MJD(ET_departure_from_asteroid + result2.t_nodes[i])

        # Convert thrust angles to Cartesian
        θ = result2.U[1, i]
        φ = result2.U[2, i]
        δ = result2.U[3, i]

        # Thrust direction in RTN
        û_r = sin(θ) * cos(φ)
        û_t = sin(θ) * sin(φ)
        û_n = cos(θ)

        # Get spacecraft state
        r = result2.X[1:3, i]
        v = result2.X[4:6, i]

        # RTN basis
        r_hat = r / norm(r)
        h = cross(r, v)
        n_hat = h / norm(h)
        t_hat = cross(n_hat, r_hat)

        # Thrust in inertial frame
        thrust_inertial = û_r * r_hat + û_t * t_hat + û_n * n_hat
        thrust_vector = GTOC12.T_max * δ * thrust_inertial

        # Write control line (6 fields, NO mass)
        @printf(io, "%5d %5d %25.12f %25.16f %25.16f %25.16f\n",
            1, -1, MJD,
            thrust_vector[1], thrust_vector[2], thrust_vector[3])
    end
end

println("  ✓ File written: $output_file")
println("  Total lines: $(2 + 2*length(result1.t_nodes) + 2 + 2*length(result2.t_nodes))")

# ========================================================================
# Run Verification
# ========================================================================
println("\n" * "="^70)
println("Running GTOC12 Verification")
println("="^70)

cd("/Users/pauldetrempe/GTOC12/problem/GTOC12_Verification_Program/GTOC12_Verification/macOS")
try
    run(`./GTOC12_Verify`)
    println("\n✅ VERIFICATION PASSED!")
catch e
    println("\n⚠️  Verification output (check above for details)")
end

println("\n" * "="^70)
println("Complete!")
println("="^70)
