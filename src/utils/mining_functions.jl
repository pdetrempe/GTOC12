function deploy_a_miner(mining_ship)
    # if a miner is being deployed, miner count decreases by 1
    # pop 1 miner off the array of miners in miners_onboard
    mining_ship.miner_count -= 1
    pop!(mining_ship.miners_onboard)
    mining_ship.mass_total = mining_ship_total_mass(mining_ship)

    return mining_ship

end


function recover_a_miner(mining_ship, asteroid_miner, current_time::Float64)
    # if a miner is being recovered, miner count increases by 1
    # add 1 miner to the array of miners in miners_onboard
    mining_ship.miner_count += 1
    append!(mining_ship.miners_onboard, asteroid_miner)
    seconds_to_years = 365 * 24 * 60 * 60
    dirt = asteroid_miner.rate_of_mining * (current_time - asteroid_miner.deployment_time) * seconds_to_years
    mining_ship.mass_collected += dirt
    mining_ship.mass_total = mining_ship_total_mass(mining_ship)

    return mining_ship

end

function mining_ship_total_mass(mining_ship)
        mass_total = sum([mining_ship.miners_onboard[i].mass for i in 1:mining_ship.miner_count]) 
            + mining_ship.mass_collected + mining_ship.mass_dry + mining_ship.mass_wet

        return mass_total

end