using GTOC12
using Plots

# Find asteroid closest (in terms of orbital energy) to the Earth
_, ID_min = findmin(abs.(GTOC12.asteroid_df.sma / au2m .- 1))
asteroid = Asteroid(ID_min)

# Try out shooting method to hit asteroid
furnish_all_kernels()
DV₀ = [-1000; 0; 0]
x₀ = get_body_state(GTOC12.Earth; ET=GTOC12.ET₀) + [0; 0; 0; DV₀[:]]

# Transfer time
Δt = 1 / 2 * 365 * 24 * 3600
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


# Try running continuous burn over this arc
# STates outputted in Cartesian 
states_out, m_out, controls_out, time_out = optimize_continuous_arc(x₀⁺, x_target, Δt; asteroid_ID=ID_min, m₀=3000.0)

# performance check 

init_state = states_out[1]
final_state = states_out[end]
pos_errori = sqrt(sum(x_target[1:3].^2)) - sqrt(sum(init_state[1:3].^2))
pos_errorf = sqrt(sum(x_target[1:3].^2)) - sqrt(sum(final_state[1:3].^2))
vel_errori = sqrt(sum(x_target[4:6].^2)) - sqrt(sum(init_state[4:6].^2))
vel_errorf = sqrt(sum(x_target[4:6].^2)) - sqrt(sum(final_state[4:6].^2))

# reshape state matrix for visualization 

states_out_matrix = mapreduce(permutedims, vcat, states_out)


plot(states_out_matrix)

plot(time_out, controls_out[1,:])


plot_body(asteroid; ETs=ETs, color=colormap("Reds"))
plot_body!(GTOC12.Earth; ETs=ETs, color=colormap("Blues"))
plot_coast!(states_out[1], Δt; label="Coast 1", color=colormap("Greens"))