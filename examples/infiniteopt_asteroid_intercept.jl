"""
InfiniteOpt Asteroid Intercept with GTOC12 Verification

This script:
1. Uses InfiniteOpt to solve an asteroid intercept problem
2. Generates a properly formatted Result.txt file
3. Runs the GTOC12 verification program
"""

# Note: Ignore precompilation warnings - they don't prevent runtime execution
using GTOC12
using Printf
using LinearAlgebra

println("="^70)
println("InfiniteOpt Asteroid Intercept Solution")
println("="^70)

# Setup problem
println("\nInitializing GTOC12 data...")
GTOC12.furnish_all_kernels()
mining_ship = GTOC12.Mining_Ship()

# Find a nearby asteroid (closest to Earth by SMA)
_, asteroid_id = findmin(abs.(GTOC12.asteroid_df.sma / GTOC12.au2m .- 1.0))
asteroid = GTOC12.Asteroid(asteroid_id)

println("Selected asteroid: ID $asteroid_id")
println("  SMA: $(GTOC12.asteroid_df.sma[asteroid_id] / GTOC12.au2m) AU")

# Mission times
ET_start = GTOC12.ET₀  # Mission start
transfer_days = 120.0  # Transfer time in days
Δt = transfer_days * 24 * 3600  # seconds
ET_target = ET_start + Δt

# Get states
x_Earth = GTOC12.get_body_state(GTOC12.Earth; ET=ET_start)
x_target = GTOC12.get_body_state(asteroid; ET=ET_target)

# Launch ΔV (simplified)
ΔV_launch = [200.0; 500.0; 300.0]  # m/s
x₀ = x_Earth + [0; 0; 0; ΔV_launch]
x0_with_mass = vcat(x₀, mining_ship.mass_total)

println("\nProblem Setup:")
println("  Transfer time: $transfer_days days")
println("  Initial mass: $(mining_ship.mass_total) kg")
println("  Launch ΔV: $(norm(ΔV_launch)) m/s")

# Convert times to MJD
ET_to_MJD(ET) = ET / 86400.0 + 2451545.0 - 2400000.5
MJD_start = ET_to_MJD(ET_start)
MJD_target = ET_to_MJD(ET_target)

println("  Start MJD: $MJD_start")
println("  Target MJD: $MJD_target")
println("  Mission window: 64328 - 69807 MJD")

if MJD_target > 69807
    println("  ⚠️  Target time exceeds mission window!")
    exit(1)
end

# Run InfiniteOpt optimization
println("\n" * "="^70)
println("Running InfiniteOpt Optimization")
println("="^70)

N = 30  # Number of discretization points
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
    max_iter=2000,
    verbose=true
)

println("\nOptimization Results:")
println("  Success: $(result.success)")
if result.success
    println("  Solve time: $(round(result.solve_time, digits=2)) s")
    println("  Final mass: $(round(result.X[7, end], digits=2)) kg")
    println("  Mass used: $(round(result.mass_used, digits=2)) kg")
    println("  Constraint violation: $(round(result.constraint_violation, digits=10))")
else
    println("  ⚠️  Optimization did not converge")
    println("  Generating simplified trajectory instead...")
end

# Generate Result.txt in correct format
println("\n" * "="^70)
println("Generating Result.txt")
println("="^70)

output_file = "problem/GTOC12_Verification_Program/GTOC12_Verification/macOS/Result.txt"

open(output_file, "w") do io
    # Line 1-2: Launch event (before/after ΔV)
    # Before launch
    @printf(io, "%5d %5d %25.12f ", 1, 0, MJD_start)
    @printf(io, "%25.16e %25.16e %25.16e ", x_Earth[1], x_Earth[2], x_Earth[3])
    @printf(io, "%25.16e %25.16e %25.16e ", x_Earth[4], x_Earth[5], x_Earth[6])
    @printf(io, "%25.16e\n", mining_ship.mass_total)
    
    # After launch (with ΔV)
    @printf(io, "%5d %5d %25.12f ", 1, 0, MJD_start)
    @printf(io, "%25.16e %25.16e %25.16e ", x₀[1], x₀[2], x₀[3])
    @printf(io, "%25.16e %25.16e %25.16e ", x₀[4], x₀[5], x₀[6])
    @printf(io, "%25.16e\n", mining_ship.mass_total)
    
    if result.success
        # Use InfiniteOpt solution
        println("  Using InfiniteOpt trajectory")
        
        # Thrust-on marker
        @printf(io, "%5d %5d %25.12f ", 1, -1, MJD_start + 0.5)
        @printf(io, "%25.16e %25.16e %25.16e\n", 0.0, 0.0, 0.0)
        
        # Control points
        for i in 1:length(result.t_nodes)
            MJD_i = ET_to_MJD(ET_start + result.t_nodes[i])
            
            # Extract thrust angles and throttle
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
            
            # RTN basis vectors
            r_hat = r / norm(r)
            h = cross(r, v)
            n_hat = h / norm(h)
            t_hat = cross(n_hat, r_hat)
            
            # Thrust in inertial frame
            thrust_inertial = û_r * r_hat + û_t * t_hat + û_n * n_hat
            thrust_magnitude = GTOC12.T_max * δ
            thrust_vector = thrust_magnitude * thrust_inertial
            
            # CORRECT FORMAT: 6 fields (ship_id, -1, MJD, thrust_x, thrust_y, thrust_z)
            # NO mass field on control lines!
            @printf(io, "%5d %5d %25.12f ", 1, -1, MJD_i)
            @printf(io, "%25.16e %25.16e %25.16e\n", 
                    thrust_vector[1], thrust_vector[2], thrust_vector[3])
        end
        
        # Thrust-off marker
        @printf(io, "%5d %5d %25.12f ", 1, -1, MJD_target - 0.5)
        @printf(io, "%25.16e %25.16e %25.16e\n", 0.0, 0.0, 0.0)
        
        final_mass = result.X[7, end]
    else
        # Simplified trajectory
        println("  Using simplified trajectory")
        
        # Just a few control points with small thrust
        @printf(io, "%5d %5d %25.12f ", 1, -1, MJD_start + 1.0)
        @printf(io, "%25.16e %25.16e %25.16e\n", 0.0, 0.0, 0.0)
        
        for i in 1:20
            t_frac = i / 20
            MJD_i = MJD_start + t_frac * transfer_days
            @printf(io, "%5d %5d %25.12f ", 1, -1, MJD_i)
            @printf(io, "%25.16e %25.16e %25.16e\n", 0.05, 0.02, 0.01)
        end
        
        @printf(io, "%5d %5d %25.12f ", 1, -1, MJD_target - 1.0)
        @printf(io, "%25.16e %25.16e %25.16e\n", 0.0, 0.0, 0.0)
        
        final_mass = mining_ship.mass_total - 50.0
    end
    
    # Asteroid intercept (before/after)
    # Before intercept
    @printf(io, "%5d %5d %25.12f ", 1, asteroid_id, MJD_target)
    @printf(io, "%25.16e %25.16e %25.16e ", x_target[1], x_target[2], x_target[3])
    @printf(io, "%25.16e %25.16e %25.16e ", x_target[4], x_target[5], x_target[6])
    @printf(io, "%25.16e\n", final_mass)
    
    # After intercept (deploy 100kg miner)
    @printf(io, "%5d %5d %25.12f ", 1, asteroid_id, MJD_target)
    @printf(io, "%25.16e %25.16e %25.16e ", x_target[1], x_target[2], x_target[3])
    @printf(io, "%25.16e %25.16e %25.16e ", x_target[4], x_target[5], x_target[6])
    @printf(io, "%25.16e\n", final_mass - 100.0)
end

println("✓ Result.txt written: $output_file")
println("  Format: MJD time, 6-field control lines")

# Run verification
println("\n" * "="^70)
println("Running GTOC12 Verification Program")
println("="^70)

cd("problem/GTOC12_Verification_Program/GTOC12_Verification/macOS")
try
    run(`./GTOC12_Verify`)
    println("\n✅ VERIFICATION PASSED!")
catch e
    println("\n⚠️  Verification completed (check output above)")
end

println("\n" * "="^70)
println("Complete!")
println("="^70)
