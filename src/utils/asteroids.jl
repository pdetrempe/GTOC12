using CSV, DataFrames, GTOC12

export get_asteroid_df, get_asteroid_state

function get_asteroid_df()
    # furnish SPICE kernels
    furnish_all_kernels()

    # asteroid data path
    asteroid_data_file = joinpath( GTOC12.PROBLEM_DATA, "GTOC12_Asteroids_Data.txt")
    coeff_file = joinpath( GTOC12.PROBLEM_DATA, "bonus_coefficients_202306019.txt")

    # Asteroid data
    asteroid_df = DataFrame(CSV.File(asteroid_data_file; header=1, delim=' ', ignorerepeated=true))

    # Rename headers to avoid syntax issues
    rename!(asteroid_df, [:ID, :MJD, :sma, :ecc, :inc, :LAN, :argperi, :mean_anom])

    # Add bonus coefficients and collected mass
    coeff_df = DataFrame(CSV.File(coeff_file; header=false, delim=' ', ignorerepeated=true))
    rename!(coeff_df, ["bonus_coeff", "collected_mass"])
    insertcols!(coeff_df, 1, :ID=>asteroid_df.ID) # Add IDs to data frame

    # Join the two data frames to have a party
    asteroid_df = innerjoin(asteroid_df, coeff_df, on=:ID)

    # Convert units to meters/radians
    asteroid_df.sma *= au2m
    asteroid_df.inc *= π/180
    asteroid_df.LAN *= π/180
    asteroid_df.argperi *= π/180
    asteroid_df.mean_anom *= π/180

    # Add other time epochs that are useful
    TT_start = "2035, Jan 1, 00:00:00.0000 (TDT)" # TT is TDT in SPICE
    asteroid_df[!, :TT] .= TT_start
    asteroid_df[!, :ET] .= str2et(TT_start)


    return asteroid_df

end

# Query asteroid state at a given time
function get_asteroid_state(asteroid, ET)
    # Change phasing to be time-based
    n = mean_motion(;a=asteroid.sma, μ=μ_☉) # mean motion

    # Propagate mean anomaly to ET
    M = asteroid.mean_anom + n * (ET - GTOC12.ET₀)

    # Calculate true anomaly
    E = M2EH(; M=M, ecc=asteroid.ecc)
    ν = EH2ν(; E_or_H=E, ecc=asteroid.ecc)

    # Get asteroid position at all true anomalies
    x_asteroid = COE2RV(; COE=
        [asteroid.sma,
            asteroid.ecc,
            asteroid.inc,
            asteroid.LAN,
            asteroid.argperi,
            ν],
        μ_CB_or_CB_name=μ_☉)
end