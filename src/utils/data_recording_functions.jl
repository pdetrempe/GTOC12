
function record_line(line_array, event::String, state, time_vector, mining_ship;
     control=nothing, rendez_flag::String="none", event_ID="ERROR")
    #new_line = DataFrame()
    #total_mass = mining_ship.mass_dry + mining_ship.mass_wet + mining_ship.mass_collected
    #line_array = []

    if event == "burn"
        # Event ID -1 
        # TODO
        event_ID = -1
        first_line = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector[1], " ") * join([0.0, 0.0, 0.0], " ") * "\n"
        line_array = vcat(line_array, first_line)

        for (time, thrust) in zip(time_vector, eachrow(control))
            line = string(mining_ship.ship_ID, " ", event_ID, " ", time, " ") * join(thrust, " ") * "\n"
            line_array = vcat(line_array, line)
        end

        last_line  = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector[end], " ") * join([0.0, 0.0, 0.0], " ") * "\n"
        line_array = vcat(line_array, last_line)


    elseif event == "rendezvous"
        # Event ID => asteroid ID 
        # Epoch, rx, ry, rz, vx, vy, vz, mass 
        # TODO  UPDATE EVENT ID
        before_line = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector, " ") * join(state, " ") * string(" ", mining_ship.mass_total) * "\n"
        line_array = vcat(line_array, before_line)

        if rendez_flag == "deploy"
            mining_ship = deploy_a_miner(mining_ship, event_ID, time_vector[end])
        else
            #mining_ship = recover_a_miner(mining_ship, event_ID, time_vector[end], asteroid_miner_dictionary)
            mining_ship = recover_a_miner(mining_ship, event_ID, time_vector[end])
        end

        after_line  = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector, " ") * join(state, " ") * string(" ", mining_ship.mass_total) * "\n"
        line_array = vcat(line_array, after_line)

    elseif event == "launch"
        # Event ID 0 
        event_ID = 0
        before_line = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector, " ") * join(state[1], " ") * string(" ", mining_ship.mass_total) * "\n"
        line_array = vcat(line_array, before_line)
        after_line  = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector, " ") * join(state[2], " ") * string(" ", mining_ship.mass_total) * "\n"
        line_array = vcat(line_array, after_line)
    elseif event == "venus_flyby"
        # Event ID -2
        event_ID = -2
        before_line = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector, " ") * join(state[1], " ") * string(" ", mining_ship.mass_total) * "\n"
        line_array = vcat(line_array, before_line)
        after_line  = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector, " ") * join(state[2], " ") * string(" ", mining_ship.mass_total) * "\n"
        line_array = vcat(line_array, after_line)
    elseif event == "earth_flyby"
        # Event ID -3
        event_ID = -3
        before_line = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector, " ") * join(state[1], " ") * string(" ", mining_ship.mass_total) * "\n"
        line_array = vcat(line_array, before_line)
        after_line  = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector, " ") * join(state[2], " ") * string(" ", mining_ship.mass_total) * "\n"
        line_array = vcat(line_array, after_line)
    elseif event == "mars_flyby"
        # Event ID -4
        event_ID = -4
        before_line = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector, " ") * join(state[1], " ") * string(" ", mining_ship.mass_total) * "\n"
        line_array = vcat(line_array, before_line)
        after_line  = string(mining_ship.ship_ID, " ", event_ID, " ", time_vector, " ") * join(state[2], " ") * string(" ", mining_ship.mass_total) * "\n"
        line_array = vcat(line_array, after_line)
    end

    return line_array, mining_ship

end