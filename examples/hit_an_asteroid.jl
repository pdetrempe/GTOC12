using GTOC12
using Plots
using TrajectoryOptimization
using RobotDynamics
using RobotDynamics: @autodiff
using StaticArrays
using ForwardDiff, FiniteDiff
using LinearAlgebra

const TO = TrajectoryOptimization
const RD = RobotDynamics

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
# TODO: add graceful error handling
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
x₀⁺, xₜ = fixed_time_single_shoot( x₀, t_transfer, r_target; print_iter=true)
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

struct Asteroid
    sma::Float64
    ecc::Float64
    inc::Float64
    LAN::Float64
    argperi::Float64
    mean_anom::Float64
    ET₀::Float64
    oe::SVector{6,Float64}
    x_ET₀::SVector{6,Float64}

    # Calculate relevant quantities at construction
    function Asteroid(sma::Float64, ecc::Float64, inc::Float64, LAN::Float64, argperi::Float64, mean_anom::Float64, ET₀::Float64)
        ν = M2ν(; M=mean_anom, ecc=ecc)
        oe = @SVector [sma, ecc, inc, LAN, argperi, ν]
        x_ET₀ = COE2RV( COE=oe )
        return new( sma, ecc, inc, LAN, argperi, mean_anom, ET₀, oe, x_ET₀)
    end
end

function Asteroid(;sma, ecc, inc, LAN, argperi, mean_anom, ET₀=GTOC12.ET₀ )
    return Asteroid( sma, ecc, inc, LAN, argperi, mean_anom, ET₀ )
end

# Initialize from ID/data frame row
function Asteroid(ID::Int)
    asteroid = GTOC12.asteroid_df[ID, :]
    return Asteroid(;sma=asteroid.sma, 
                    ecc=asteroid.ecc, 
                    inc=asteroid.inc, 
                    LAN=asteroid.LAN, 
                    argperi=asteroid.argperi, 
                    mean_anom=asteroid.mean_anom, 
                    ET₀=asteroid.ET )
end



@autodiff struct ImpulsiveSpacecraft <: TO.DiscreteDynamics
    ET_launch::Float64
    target_asteroid::Asteroid
end

function ImpulsiveSpacecraft(;ET_launch, target_asteroid)
    return ImpulsiveSpacecraft(ET_launch, target_asteroid)
end

function Base.copy(c::ImpulsiveSpacecraft)
    ImpulsiveSpacecraft(c.ET_launch, c.target_asteroid)
end

function RD.discrete_dynamics(model::ImpulsiveSpacecraft, x, u, t, dt)
    ΔV = u
    x₀⁺ = x + [0;0;0; ΔV[:] ]
    xₜ =  propagate_universal(x₀⁺, dt)
    SVector{length(xₜ)}(xₜ)
end

RD.state_dim(::ImpulsiveSpacecraft) = 6
RD.control_dim(::ImpulsiveSpacecraft) = 3

## Problem definition

# Model and discretization
target_asteroid = Asteroid(ID_min)
model = ImpulsiveSpacecraft(GTOC12.ET₀, target_asteroid)
n,m = RD.dims(model)
tf = t_transfer  # final time (sec)
N = 3           # number of knot points, 2 for simple start/end impulses
dt = tf / (N-1)  # time step (sec)
# @test (n,m) == (6,3)

# Objective
# Note: May have to tune on a segment-by-segment basis
x0 = SVector{length(x₀)}(x₀)# initial state
xf = SVector{length(x_target)}(x_target) # final state

# TODO: Convert costs to just Σ(ΔV)
Q = Diagonal(@SVector ones(n))
R = Diagonal(@SVector ones(m))
Qf = Diagonal( @SVector ones(n))
obj = LQRObjective(Q, R, Qf, xf, N)

# Add constraints
conSet = ConstraintList(n,m,N)
add_constraint!(conSet, GoalConstraint(xf), N)
bnd = BoundConstraint(n,m, u_min=-6000/sqrt(m)*ones(m), u_max=6000/sqrt(m)*ones(m)) # 6 km/s limit on hyperbolic excess velocity
add_constraint!(conSet, bnd, 1:N-1)

# Initial guess
# u0 = SVector{length(ΔV_departure_1)}(ΔV_departure_1) # Try warm-starting with shooting method solution?
u0 = @SVector zeros(m)

# Set up problem

# TODO: Change problem formulation to constrain final position and add final velocity difference to COST

prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)
initial_controls!(prob, u0)
rollout!(prob);


# Uh, solve the problem?
using Altro
opts = SolverOptions(;
    constraint_tolerance = 1.0,
    cost_tolerance = 1.0,
    cost_tolerance_intermediate = 1.0,
    max_cost_value = 1.0e30,
    max_state_value = 1.0e30,
    max_control_value = 6.0e30
    # Optimality Tolerances
    

    # iLQR
    # expected_decrease_tolerance = 1e-10
    # iterations_inner = 300
    # dJ_counter_limit = 10
    # square_root = false
    # line_search_lower_bound = 1e-8
    # line_search_upper_bound = 10.0
    # line_search_decrease_factor = 0.5
    # iterations_linesearch = 20

    # static_bp = true
	# save_S = false
    # closed_loop_initial_rollout = false
)

solver = ALTROSolver(prob, opts);
solve!(solver)
println("Cost: ", cost(solver))
# println("Constraint violation: ", max_violation(solver))
println("Iterations: ", iterations(solver))

# Extract states and controls
X = states(prob)
U = controls(prob)
t = gettimes(prob)

# Xrollout = [copy(x0) for k = 1:N]
# Urollout = [copy(u0) for k = 1:N]
# for k = 1:N-1
#     Xrollout[k+1] = RD.discrete_dynamics(
#         get_model(prob, k), Xrollout[k], Urollout[k], dt*(k-1), dt
#     )
# end
# # @test Xrollout ≈ X 


# # Convert to matrices
# Xmat = hcat(Vector.(X)...)
# Umat = hcat(Vector.(U)...)
