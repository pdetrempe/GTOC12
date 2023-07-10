using GTOC12
using Plots
using LinearAlgebra
using OrdinaryDiffEq
using ForwardDiff

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1.))
asteroid = Asteroid(ID_min)

# Try out shooting method to hit asteroid
furnish_all_kernels()
DV₀ = [-1000.; 0.; 3500.]
x₀ = get_body_state(GTOC12.Earth; ET=GTOC12.ET₀) + [0; 0; 0; DV₀[:]]

# Transfer time
Δt = 0.75 * 365 * 24 * 3600
ΔV∞_max = 6000.0;

# Fixed-time shooting to hit asteroid
# Get problem parameters: asteroid position (target), initial state, transfer time
ET_target = GTOC12.ET₀ + Δt
x_target = get_body_state(asteroid; ET=ET_target)
r_target = view(x_target, 1:3)

# NOTE: This is super brittle to initial conditions
# TODO: add a Lambert targetter for doing initial guess
# TODO: do a time sweep of Lambert trajectories to minimize DV
x₀⁺, xₜ = fixed_time_single_shoot(x₀, Δt, r_target; print_iter=false)
ΔV_departure_1 = (x₀⁺-x₀)[4:6]
ΔV_arrival_1 = (xₜ-x_target)[4:6]

# Plot Earth/asteroid/spacecraft
ETs = [GTOC12.ET₀, ET_target]
# plot_asteroid_from_df_row(asteroid; ET_in=ETs)
plot_body(asteroid; ETs=ETs, color=colormap("Reds"))
plot_body!(GTOC12.Earth; ETs=ETs, color=colormap("Blues"))
plot_coast!(x₀⁺, Δt; label="Coast 1", color=colormap("Greens"))

# TODO: Optimize the above for DV/via Time-of-flight

# TODO: For launch case, just need initial time, not states
#states_out, controls_out, time_out = optimize_impulsive_launch(x₀, x_target, Δt; ΔV₀=ΔV_departure_1, asteroid_ID=ID_min)

# Delta V from impulsive 
  # 1.1632731732968505e-6
  # 2.408266604131654e-6
  #-4.037442985551002e-6
## Mining Ship Mass Constraint
# m0 = md + mp + I*ms <= 3000 kg
# m0 = 500kg + mp + 20*40 = 1300kg + mp
# md = 500 kg, dry mass 
# mp = prop mass
# ms = miner mass, 40 kg
# I = # miners <= 20

#----------------- Try Boundary condition version
function getAB(x)
    p = x[1];
    f = x[2];
    g = x[3];
    h = x[4];
    k = x[5];
    l = x[6]
    q = 1.0 + f*cos(l) + g*sin(l)
    s = sqrt(1.0 + h^2 + k^2)
    A = [0, 0, 0, 0, 0, sqrt(μ*p*(q/p)^2)]

    B = [0                2*p/q*sqrt(p/μ)                0;
         sqrt(p/μ)*sin(l) sqrt(p/μ)*1.0/q*(q+1)*cos(l)+f  -sqrt(p/μ)*g/q*(h*sin(l)-k*cos(l));
        -sqrt(p/μ)*cos(l) sqrt(p/μ)*1.0/q*(q+1)*sin(l)+g   sqrt(p/μ)*g/q*(h*sin(l)-k*cos(l));
         0                0                              sqrt(p/μ)*s*cos(l)/(2*q);
         0                0                              sqrt(p/μ)*s*sin(l)/(2*q);
         0                0                              sqrt(p/μ)*1.0/q*(sin(l)-k*cos(l))]

    return A, B
end

function low_thrust_optimal_control!(dstate, state, p, t)
    # unpack parameters
    # no parameters to unpack at this time, could unpack a_grav here (and T i suppose)?
    #a_grav, T = p
    T = 0.6 # Max Thrust 
    a_grav = [0.,0.,0.]

    # unpack state
    x = state[1:6]; m = state[7]; λ = state[8:13]

    # Define A and B matrices (time varying)
    # *************************************************
    A, B = getAB(x)

    # Optimal Control Strategy 
    # *************************************************
    #u_star = [0,0,0]
    #δ_star = 0
    u_star = -B'*λ/norm(B'*λ)  
    S = norm(B'*λ) .- 1.0        
    if S > 0
        δ_star = 1.0
    else
        δ_star = 0.0
    end

    # Spacecraft Dynamics 
    # *************************************************
    dx = A + B*a_grav + T*δ_star/m * B*u_star

    # Mass variation
    # *************************************************
    # exhaust velocity
    g0 = 9.81    # m/s^2
    Isp = 4000.0   # seconds 
    c = Isp * g0 # specific impulse and grav accel at sea level (m/s)
    #c = 1 
    dm = -T/c * δ_star

    # Costate diff eqs
    param = m, λ, δ_star, u_star, a_grav, T 
    dH_dx = ForwardDiff.gradient(x -> calculate_hamiltonian(x, param ), x) #[1:6])
    dλ = -dH_dx'

    dstate = [dx, dm, dλ]

end

function calculate_hamiltonian(x, p) 
    m, λ, δ_star, u_star, a_grav, T = p
    A, B = getAB(x)
    H = T/m*δ_star +    λ'*A    +   λ'*B*(T/m*δ_star*u_star + a_grav)
    return H
end



function bc2!(residual, state, p, t) 
    # u[1] is the beginning of the time span, and u[end] is the ending
    # TODO update params to input p_0, p_f
    x0, xf = p
    MEE_init = Cartesian2MEE(x0; μ=GTOC12.μ_☉)
    MEE_target = Cartesian2MEE(xf; μ=GTOC12.μ_☉)

    p_0 = MEE_init
    p_f = MEE_target
    # initial boundary value 
    # ***********************************
    residual[1] = state[1][1] - p_0[1] 
    residual[2] = state[1][2] - p_0[2] 
    residual[3] = state[1][3] - p_0[3] 
    residual[4] = state[1][4] - p_0[4] 
    residual[5] = state[1][5] - p_0[5] 
    residual[6] = state[1][6] - p_0[6] 
    # final boundary value 
    # ***********************************
    residual[7]  = state[end][1] - p_f[1] 
    residual[8]  = state[end][2] - p_f[2] 
    residual[9]  = state[end][3] - p_f[3] 
    residual[10] = state[end][4] - p_f[4] 
    residual[11] = state[end][5] - p_f[5] 
    residual[12] = state[end][6] - p_f[6] 
end

m = 500.0 # kg
μ=GTOC12.μ_☉
MEE_init = Cartesian2MEE(x₀⁺; μ=GTOC12.μ_☉)
state_init = vcat(MEE_init, m, zeros(6))
tspan = (0.0,Δt)
p = (x₀⁺,x_target)
bvp2 = TwoPointBVProblem(low_thrust_optimal_control!, bc2!, state_init, tspan, p)
sol2 = solve(bvp2, Vern7()) # we need to use the MIRK4 solver for TwoPointBVProblem
plot(sol2.u)


# # Try running continuous burn over this arc
# # STates outputted in Cartesian 
# # states_out, m_out, controls_out, time_out = optimize_continuous_arc(x₀⁺, x_target, Δt; asteroid_ID=ID_min, m₀=3000.0)

# # performance check 

# init_state = states_out[1]
# final_state = states_out[end]
# pos_errori = norm(x_target[1:3] - init_state[1:3])
# pos_errorf = norm(x_target[1:3] - final_state[1:3])
# vel_errori = norm(x_target[4:6] - init_state[4:6])
# vel_errorf = norm(x_target[4:6] - final_state[4:6])

# # reshape state matrix for visualization 

# states_out_matrix = mapreduce(permutedims, vcat, states_out)


# plot(states_out_matrix)

# plot(time_out, controls_out[1,:])


# plot_body(asteroid; ETs=ETs, color=colormap("Reds"))
# plot_body!(GTOC12.Earth; ETs=ETs, color=colormap("Blues"))
# plot_coast!(states_out[1], Δt; label="Coast 1", color=colormap("Greens"))