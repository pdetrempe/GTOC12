function deploy_a_miner(mining_ship, event_ID, deployment_time)
    # if a miner is being deployed, miner count decreases by 1
    # pop 1 miner off the array of miners in miners_onboard
    mining_ship.miner_count -= 1
    #pop!(mining_ship.miners_onboard)
    mining_ship.mass_total = mining_ship_total_mass(mining_ship)
    if haskey(mining_ship.asteroid_miner_dictionary,"$event_ID")
        mining_ship.asteroid_miner_dictionary["$event_ID"] = [mining_ship.asteroid_miner_dictionary["$event_ID"], deployment_time]
    else
        mining_ship.asteroid_miner_dictionary["$event_ID"] = deployment_time
    end

    return mining_ship

end


function recover_a_miner(mining_ship, event_ID, recovery_time::Float64)
    # if a miner is being recovered, miner count increases by 1
    # add 1 miner to the array of miners in miners_onboard
    # for loop to protect for edge case with multiple miners on 1 asteroid 
    seconds_to_years = 365 * 24 * 60 * 60
    for deployed_time in mining_ship.asteroid_miner_dictionary["$event_ID"]

        mining_ship.miner_count += 1
    
        dirt = mining_ship.rate_of_mining * (recovery_time - deployed_time) * seconds_to_years
        mining_ship.mass_collected += dirt
    end
    mining_ship.mass_total = mining_ship_total_mass(mining_ship)

    return mining_ship

end

function mining_ship_total_mass(mining_ship)
        mass_total = mining_ship.miner_mass*mining_ship.miner_count + mining_ship.mass_collected + mining_ship.mass_dry + mining_ship.mass_wet
        return mass_total
end