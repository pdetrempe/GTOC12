
mutable struct asteroid_miner
    rate_of_mining::Float64                       # rate of mining kg/year
    mass::Float64                                 # mass of miner, kg
    deployment_time::Float64                      # time when deployed, sec

    function asteroid_miner(;
        rate_of_mining = 10.0,                    # (kg/year)
        mass = 40.0,                              # (kg)
        deployment_time = 0.0,                    # ET
        )
        new(rate_of_mining, mass, deployment_time)
    end
end

mutable struct Mining_Ship
    ship_ID::Int
    miner_count::Float64 
    miners_onboard::Array                         # array of miners on board
    mass_dry::Float64                             # Dry mass (kg)
    mass_wet::Float64                             # wet mass from prop (kg)
    mass_collected::Float64                       # total collected mass on board Mining Ship
    mass_total::Float64                           # total mass of Mining Ship 
    current_state::Array                          # current orbital state, Cartesian Coordinates

    function Mining_Ship(;
        ship_ID = 1,
        miner_count = 20,                                   # (cnt)
        miners_onboard = [GTOC12.asteroid_miner() for i in 1:miner_count],  # (-)
        mass_dry = 500.,                                    # (kg) 
        mass_wet = 1300.,                                   # (kg)
        mass_collected = 0.0,                               # (kg)
        mass_total = sum([miners_onboard[i].mass for i in 1:miner_count]) + mass_collected + mass_dry + mass_wet,  # (kg)
        current_state = zeros(6),                           # [r, v]
        )
        new(ship_ID, miner_count, miners_onboard, mass_dry, mass_wet, mass_collected, mass_total, current_state)
    end
end