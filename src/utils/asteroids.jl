using CSV, DataFrames, GTOC12

export get_asteroid_df

function get_asteroid_df()
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

end