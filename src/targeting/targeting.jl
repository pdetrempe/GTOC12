using ForwardDiff

export fixed_time_single_shoot, calculate_miss_distance

# Get STM: use ForwardDiff to get dr_f/dv_0
# Actual "cost" is miss distance
# Need f(v₀) = dr to fit ForwardDiff format
function calculate_miss_distance(x_aug)
    # Initial guess state values
    x₀ = x_aug[1:6]
    t = x_aug[7]
    rₜ_des = x_aug[8:10]

    xₜ = propagate_universal(x₀, t)

    cost = xₜ[1:3] - rₜ_des
end

# TODO: Think about clever way to package problem into a type
# such that can have one giant vector for calculating adjoint_sensitivities

function fixed_time_single_shoot( x₀, t_shoot, r_target; MAX_ITER = 50, print_iter=false)

    num_iter = 0
    x₀⁺ = copy(x₀)
    drₜ = x₀[1:3] # Initialize miss distance
        while num_iter < MAX_ITER && norm(drₜ) > GTOC12.POS_ABS_TOL
            # Form vector used by ForwardDiff problem
            x_aug = vcat(x₀⁺, t_shoot, r_target)

            # Calculate miss distance from current propagation
            drₜ = calculate_miss_distance(x_aug)
            if print_iter
                println(drₜ)
            end

            # Calculate sensitivity of solution to all parameters
            drₜ_dx_aug = ForwardDiff.jacobian(calculate_miss_distance, x_aug)

            # Get partials w.r.t. design variables
            # (In this case, initial velocity and time of flight)
            ∂rₜ_∂v₀ = drₜ_dx_aug[:, 4:6] # Variation of final position w.r.t. initial velocity
            # ∂rₜ_∂TOF = dr_dx_aug[:, 7]  # Variation of final position w.r.t. time-of-flight

            # ^^^TODO: handle indices more intelligently for packing/unpacking


            # Update initial state via Newton's method
            # See Pavlak Thesis Eq. 3.12
            dv₀ = -∂rₜ_∂v₀\drₜ
            # dTOF = -∂rₜ_∂TOF\drₜ

            x₀⁺[4:6] += dv₀

            # Variable time stuff isn't working interestingly
            # global t_shoot += dTOF

            num_iter += 1

        end
    xₜ = propagate_universal(x₀⁺, t_shoot)
    return x₀⁺, xₜ
end
