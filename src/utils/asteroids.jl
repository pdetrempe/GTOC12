using CSV, DataFrames, SPICE, GTOC12

export get_asteroid_df, get_asteroid_state, asteroid_df

function get_asteroid_df()
    # furnish SPICE kernels if they haven't been loaded yet
    furnish_all_kernels()

    # asteroid data path
    asteroid_data_file = joinpath( PROBLEM_DATA, "GTOC12_Asteroids_Data.txt")
    coeff_file = joinpath( PROBLEM_DATA, "bonus_coefficients_202306019.txt")

    # Asteroid data
    asteroid_df = CSV.read(asteroid_data_file, DataFrame; header=1, delim=' ', ignorerepeated=true)

    # Rename headers to avoid syntax issues
    rename!(asteroid_df, [:ID, :MJD, :sma, :ecc, :inc, :LAN, :argperi, :mean_anom])

    # Add bonus coefficients and collected mass
    coeff_df = CSV.read(coeff_file, DataFrame; header=false, delim=' ', ignorerepeated=true)
    rename!(coeff_df, ["bonus_coeff", "collected_mass"])

    # Join the two data frames to have a party
    asteroid_df = hcat(asteroid_df,coeff_df) # We assume the dataframes are in the same order by ID

    # Convert units to meters/radians
    asteroid_df.sma *= au2m
    asteroid_df.inc = map(deg2rad, asteroid_df.inc)
    asteroid_df.LAN = map(deg2rad, asteroid_df.LAN)
    asteroid_df.argperi = map(deg2rad, asteroid_df.argperi)
    asteroid_df.mean_anom = map(deg2rad, asteroid_df.mean_anom)

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

# Only load/manipulate dataframe once
const asteroid_df = get_asteroid_df()

# Problem start time
const ET₀ = asteroid_df[1, :ET]