function deploy_a_miner(mining_ship::Mining_Ship)
    # if a miner is being deployed, miner count decreases by 1
    # pop 1 miner off the array of miners in miners_onboard
    mining_ship.miner_count -= 1
    pop!(mining_ship.miners_onboard)

    return mining_ship

end


function recover_a_miner(mining_ship::Mining_Ship, asteroid_miner::miner, current_time::Float64)
    # if a miner is being recovered, miner count increases by 1
    # add 1 miner to the array of miners in miners_onboard
    mining_ship.miner_count += 1
    append!(mining_ship.miners_onboard, asteroid_miner)
    seconds_to_years = 365 * 24 * 60 * 60
    dirt = asteroid_miner.rate_of_mining * (current_time - asteroid_miner.deployment_time) * seconds_to_years
    mining_ship.mass_collected += dirt

    return mining_ship

end

