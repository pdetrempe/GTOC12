using GTOC12
using SPICE
using ForwardDiff
using LinearAlgebra
using Plots
using DifferentialEquations

function getAB(x)
    p = x[1];
    f = x[2];
    g = x[3];
    h = x[4];
    k = x[5];
    l = x[6]
    q = 1 + f*cos(l) + g*sin(l)
    s = sqrt(1 + h^2 + k^2)
    A = [0, 0, 0, 0, 0, sqrt(μ*p*(q/p)^2)]

    B = [0                2*p/q*sqrt(p/μ)                0;
         sqrt(p/μ)*sin(l) sqrt(p/μ)*1/q*(q+1)*cos(l)+f  -sqrt(p/μ)*g/q*(h*sin(l)-k*cos(l));
        -sqrt(p/μ)*cos(l) sqrt(p/μ)*1/q*(q+1)*sin(l)+g   sqrt(p/μ)*g/q*(h*sin(l)-k*cos(l));
         0                0                              sqrt(p/μ)*s*cos(l)/(2*q);
         0                0                              sqrt(p/μ)*s*sin(l)/(2*q);
         0                0                              sqrt(p/μ)*1/q*(sin(l)-k*cos(l))]

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
    #xₜ = propagate_keplerian(x₀, t_shoot)

    # Optimal Control Strategy 
    # *************************************************
    #u_star = [0,0,0]
    #δ_star = 0
    # TODO add logic to make this not a discontinuous change
    u_star = -B'*λ/norm(B'*λ)  
    S = norm(B'*λ) .- 1        
    if S > 0
        δ_star = 1
    else
        δ_star = 0
    end

    # Spacecraft Dynamics 
    # *************************************************
    dx = A + T*δ_star/m * B*u_star  #B*a_grav 

    # Mass variation
    # *************************************************
    # exhaust velocity
    g0 = 9.81    # m/s^2
    Isp = 4000   # seconds 
    c = Isp * g0 # specific impulse and grav accel at sea level (m/s)
    #c = 1 
    dm = -T/c * δ_star

    # Costate diff eqs
    # *************************************************
    p = m, λ, δ_star, u_star, a_grav, T 
    dH_dx = ForwardDiff.gradient(x -> calculate_hamiltonian(x, p ), x) #[1:6])
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
    p_0, p_f = p
    #cart_init = get_planet_state("EARTH", GTOC12.ET₀) # + [0;0;0;DV₀[:]]

    #cart_target = [-1.142626400391252e11, -7.363183565567313e10, 3.999306363734976e10, 2629.801421941236, -24735.59255707686, -2945.5498712683434]
    #MEE_init = Cartesian2MEE(cart_init, μ=GTOC12.μ_☉)
    #cart_target = [-2.526739015640893e10, 1.4491856083257202e11, -1.1774074660398066e7, -35830.39920517383, -5218.801143532724, -0.3853660712700435]
    MEE_0 = Cartesian2MEE(p_0, μ=GTOC12.μ_☉)
    MEE_f = Cartesian2MEE(p_f, μ=GTOC12.μ_☉)
    #p_0 = MEE_init
    #p_f = MEE_target
    # initial boundary value 
    # ***********************************
    residual[1] = state[1][1] - MEE_0[1] 
    residual[2] = state[1][2] - MEE_0[2] 
    residual[3] = state[1][3] - MEE_0[3] 
    residual[4] = state[1][4] - MEE_0[4] 
    residual[5] = state[1][5] - MEE_0[5] 
    residual[6] = state[1][6] - MEE_0[6] 
    # final boundary value 
    # ***********************************
    residual[7]  = state[end][1] - MEE_f[1] 
    residual[8]  = state[end][2] - MEE_f[2] 
    residual[9]  = state[end][3] - MEE_f[3] 
    residual[10] = state[end][4] - MEE_f[4] 
    residual[11] = state[end][5] - MEE_f[5] 
    residual[12] = state[end][6] - MEE_f[6] 
end


# Import asteroids
asteroid_df = get_asteroid_df()

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(asteroid_df.sma .- 1))

# Plot the asteroid trajectory
#asteroid = asteroid_df[ID_min, :]
#plot_asteroid_from_df_row(asteroid)

# Plot Earth
#plot_planet!(planet="EARTH"; label="Earth", color=colormap("Blues"))

# Try out shooting method to hit asteroid

DV₀ = [-GTOC12.v∞_max; 0; 0]
#kep_init = get_planet_state("EARTH", GTOC12.ET₀) + [0;0;0;DV₀[:]]
#MEE_init = keplerian2MEE(a=kep_init[1], 
#                         e=kep_init[2], 
#                         i=kep_init[3], 
#                         Ω=kep_init[4], 
#                         ω=kep_init[5], 
#                         ν=kep_init[6])
cart_init = get_planet_state("EARTH", GTOC12.ET₀) #  + [0;0;0;DV₀[:]]
MEE_init  = Cartesian2MEE(cart_init, μ=GTOC12.μ_☉)
m = 500 # kg
μ=GTOC12.μ_☉
state_init = vcat(MEE_init, m, zeros(6))
t_shoot = 1 / 2 * 365 * 24 * 3600
#t_shoot = 3600
#t_shoot = 30
t_shoot = 1 / 2 * 365 * 24 * 3600
tspan = (0, t_shoot)

#sol2 = solve(bvp2, MIRK4(), dt=.1) # we need to use the MIRK4 solver for TwoPointBVProblem
#plot(sol2)

#state = state_init
#x = state[1:6]; m = state[7]; λ = state[8:13]
#δ_star = 0; u_star = [0,0,0]; T = 3 # N
#a_grav = u_star
#p = (m, λ, δ_star, u_star, a_grav, T)
#H = calculate_hamiltonian(x, p) # m, λ, δ_star, u_star, a_grav, T )
#dH_dx = ForwardDiff.jacobian(x => calculate_hamiltonian(x, m, λ, δ_star, u_star, a, T ), x)
#dH_dx = ForwardDiff.jacobian(x -> calculate_hamiltonian(x, p ), x[1:6])
#dH_dx = ForwardDiff.gradient(x -> calculate_hamiltonian(x, p ), x) #[1:6])
#calculate_hamiltonian( m, λ, δ_star, u_star, a_grav, T, x)


#cart_init = get_planet_state("EARTH", GTOC12.ET₀) + [0;0;0;DV₀[:]]
#function Cartesian2MEE(x⃗; μ)
#MEE_init = Cartesian2MEE(cart_init, μ=GTOC12.μ_☉)


t_shoot = Float64(1 / 2 * 365 * 24 * 3600)
tspan = (0., t_shoot)
x0 = get_planet_state("EARTH", GTOC12.ET₀) #  + [0;0;0;DV₀[:]]
xf = propagate_universal(x0, t_shoot) #; μ=GTOC12.μ_☉, tol=1e-6)
p = (x0, xf)
bvp2 = TwoPointBVProblem(low_thrust_optimal_control!, bc2!, state_init, tspan, p)

sol2 = solve(bvp2) #, MIRK4(), dt=.1) # we need to use the MIRK4 solver for TwoPointBVProblem
