
function record_line(event::String, state0, statef, current_time , mining_ship, control; rendez_flag::String="none")
    #new_line = DataFrame()
    #total_mass = mining_ship.mass_dry + mining_ship.mass_wet + mining_ship.mass_collected
    line_array = []

    if event == "burn"
        # Event ID -1 
        # TODO


    elseif event == "rendezvous"
        # Event ID => asteroid ID 
        # Epoch, rx, ry, rz, vx, vy, vz, mass 
        before_line = string(mining_ship.ship_ID, " ", event_ID, " ") * join(state, " ") * mining_ship.mass_total
        line_array.append!(before_line)

        if rendez_flag == "deploy"
            mining_ship = deploy_a_miner(mining_ship)
        else
            mining_ship = recover_a_miner(mining_ship, miner, current_time)
        end

        after_line  = string(mining_ship.ship_ID, " ", event_ID, " ") * join(state, " ") * mining_ship.mass_total
        line_array.append!(after_line)

    elseif event == "launch"
        # Event ID 0 
        event_ID = 0
        before_line = string(mining_ship.ship_ID, " ", event_ID, " ") * join(state0, " ") * string(mining_ship.mass_total)
        append!(line_array, before_line)
        after_line  = string(mining_ship.ship_ID, " ", event_ID, " ") * join(statef, " ") * string(mining_ship.mass_total)
        append!(line_array, after_line)
    elseif event == "venus_flyby"
        # Event ID -2
        event_ID = -2
        before_line = string(mining_ship.ship_ID, " ", event_ID, " ") * join(state0, " ") * string(mining_ship.mass_total)
        append!(line_array, before_line)
        after_line  = string(mining_ship.ship_ID, " ", event_ID, " ") * join(statef, " ") * string(mining_ship.mass_total)
        append!(line_array, after_line)
    elseif event == "earth_flyby"
        # Event ID -3
        event_ID = -3
        before_line = string(mining_ship.ship_ID, " ", event_ID, " ") * join(state0, " ") * string(mining_ship.mass_total)
        append!(line_array, before_line)
        after_line  = string(mining_ship.ship_ID, " ", event_ID, " ") * join(statef, " ") * string(mining_ship.mass_total)
        append!(line_array, after_line)
    elseif event == "mars_flyby"
        # Event ID -4
        event_ID = -4
        before_line = string(mining_ship.ship_ID, " ", event_ID, " ") * join(state0, " ") * string(mining_ship.mass_total)
        append!(line_array, before_line)
        after_line  = string(mining_ship.ship_ID, " ", event_ID, " ") * join(statef, " ") * string(mining_ship.mass_total)
        line_array.append!(after_line)
        append!(line_array, after_line)

    end

    return line_array

end