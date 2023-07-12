function record_line(df::DataFrame, event::String, current_state, control, Mining_Ship)
    new_line = DataFrame()
    total_mass = Mining_Ship.mass_dry + Mining_Ship.mass_wet + Mining_Ship.mass_collected

    if event == "burn"
        # Event ID -1 


    elseif event == "rendezvous"
        # Event ID => asteroid ID 
        # Epoch, rx, ry, rz, vx, vy, vz, mass 
        before_line = string($Mining_Ship.ship_ID, " ", event_ID, " ") * join(state, " ") * total_mass

        after_line  = string($Mining_Ship.ship_ID, " ", event_ID, " ") * join(state, " ") * total_mass

    elseif event == "launch"
        # Event ID 0 
        event_ID = 0
        before_line = string($Mining_Ship.ship_ID, " ", event_ID, " ") * join(state, " ") * total_mass
        after_line  = string($Mining_Ship.ship_ID, " ", event_ID, " ") * join(state, " ") * total_mass
    elseif event == "venus_flyby"
        # Event ID -2
        event_ID = -2
        new_line = string($Mining_Ship.ship_ID, " ", event_ID, " ") * join(state, " ") * total_mass
    elseif event == "earth_flyby"
        # Event ID -3
        event_ID = -3
        new_line = string($Mining_Ship.ship_ID, " ", event_ID, " ") * join(state, " ") * total_mass
    elseif event == "mars_flyby"
        # Event ID -4
        event_ID = -4
        new_line = string($Mining_Ship.ship_ID, " ", event_ID, " ") * join(state, " ") * total_mass

    end




end