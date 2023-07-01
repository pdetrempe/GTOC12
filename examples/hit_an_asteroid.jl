using GTOC12
using Plots
using StaticArrays
using TrajectoryOptimization
using ForwardDiff, FiniteDiff
using LinearAlgebra

const TO = TrajectoryOptimization

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma/au2m .- 1))
asteroid = GTOC12.asteroid_df[ID_min, :]

# Try out shooting method to hit asteroid
furnish_all_kernels()
DV₀ = [-1000; 0; 0]
x₀ = get_planet_state("EARTH", GTOC12.ET₀) + [0;0;0;DV₀[:]]

# Transfer time
t_transfer = 1/2 * 365 * 24 * 3600

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = GTOC12.ET₀ + t_transfer
x_target = get_asteroid_state(asteroid, ET_target)
r_target = view(x_target, 1:3)

# NOTE: This is super brittle to initial conditions
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
x₀⁺, xₜ = fixed_time_single_shoot( x₀, t_transfer, r_target; print_iter=false)
ΔV_departure_1 = ( x₀⁺ - x₀ )[4:6]
ΔV_arrival_1 = (xₜ - x_target )[4:6]

# Add small coast period
t_coast = 24*3600

# # Try to hit Earth again
# x₀ = xₜ
# t_transfer2 = t_transfer

# ET_arrival = GTOC12.ET₀ + t_transfer + t_coast + t_transfer2
# x_target = get_planet_state("EARTH", ET_arrival)
# r_target = x_target[1:3]
# x₀⁺_2, xₜ = fixed_time_single_shoot( x₀, t_transfer, r_target; print_iter=false)
# ΔV_departure_2 = ( x₀⁺_2 - x₀ )[4:6]
# ΔV_arrival_2 = (xₜ - x_target )[4:6]


# Plot Earth/asteroid/spacecraft
ETs = [GTOC12.ET₀, ET_target]
plot_asteroid_from_df_row(asteroid; ET_in=ETs)
plot_planet!(planet="EARTH"; ET_in=ETs, label="Earth", color=colormap("Blues"))
plot_coast!(x₀⁺, t_transfer; label="Coast 1", color=colormap("Greens"))
# plot_coast!(x₀⁺_2, t_transfer2; label="Coast 2", color=colormap("Greens"))

# TODO: Optimize the above for DV

# Holy shit, this is actually a feasible solution we can use for scoring
# Let's try to write it to a scorable file

# TODO: Create an abstract type for both planets and asteroids
# Since they're both on Keplerian rails, but planets can conduct flybys and asteroids can provide resources


## Problem definition
# roblem scaling via canonical units (this makes the solver not hate you)
r⃗₀ = x₀[1:3]
v⃗₀ = x₀[4:6]
CDU = norm(r⃗₀)  # Canonical Distance Unit
CTU = √(CDU^3/GTOC12.μ_☉)# Canonical Time Unit
μ_canonical = 1.0

r⃗₀ = r⃗₀/CDU
v⃗₀ = v⃗₀/(CDU/CTU)
x₀ = vcat(r⃗₀, v⃗₀)
x_target = vcat( x_target[1:3]/CDU, x_target[4:6]/(CDU/CTU))
ΔV∞_max = 6000/(CDU/CTU) # Max of 6 km/s hyperbolic excess
ΔV_departure_1 /= (CDU/CTU)
Δt = t_transfer/CTU


# Model and discretization
target_asteroid = Asteroid(ID_min)
model = ImpulsiveSpacecraft(ET_launch=GTOC12.ET₀, target_asteroid=target_asteroid, μ=μ_canonical)
n,m = RD.dims(model)
tf = Δt  # final time (sec)
N = 2          # number of knot points, 2 for simple start/end impulses
dt = tf / (N-1)  # time step (sec)

# Objective
# Note: May have to tune on a segment-by-segment basis
x0 = SVector{length(x₀)}(x₀)# initial state
xf = SVector{length(x_target)}(x_target) # final state

# TODO: Convert costs to just Σ(ΔV)
Q = Diagonal(@SVector zeros(n)) # Don't penalize intermediate states at all.
R = 1e-10 * Diagonal(@SVector ones(m))  # Don't Penalize ΔV (control) since we get it for free from launch vehicle
Qf = Diagonal( @SVector [0;0;0;1;1;1]) # Only penalize final velocity (treat ΔV to match at the end as part of objective)
obj = LQRObjective(Q, R, Qf, xf, N)

# Add constraints
conSet = ConstraintList(n,m,N)
add_constraint!(conSet, GoalConstraint(xf, 1:3), N) # Just constrain the position at the end
bnd = BoundConstraint(n,m, u_min=-ΔV∞_max/sqrt(m)*ones(m), u_max=ΔV∞_max/sqrt(m)*ones(m)) # 6 km/s limit on hyperbolic excess velocity
add_constraint!(conSet, bnd, 1:N-1)

# Initial guess
# u0 = SVector{length(ΔV_departure_1)}(ΔV_departure_1) # Try warm-starting with shooting method solution?
u0 = @SVector zeros(m)

# Set up problem
prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)
initial_controls!(prob, u0)
rollout!(RD.InPlace(), prob);


# Uh, solve the problem?
using Altro
# See options here: https://github.com/RoboticExplorationLab/Altro.jl#list-of-options
opts = SolverOptions(;
    penalty_initial = 100.0,
    penalty_scaling = 1.0,
    constraint_tolerance = 1e-6,
    cost_tolerance = 1e-6
)

solver = ALTROSolver(prob, opts);
solve!(solver)

# Extract states and controls
X_canonical = states(prob)
U_canonical = controls(prob)
t_canonical = gettimes(prob)

# Redimensionalize output based on canonical units
X_dimensional = [vcat(X[1:3]*CDU, X[4:6]*CDU/CTU) for X in X_canonical]
U_dimensional = [U*CDU/CTU for U in U_canonical]

