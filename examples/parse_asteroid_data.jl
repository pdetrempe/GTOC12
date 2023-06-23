# using GTOC12
using CSV, DataFrames

asteroid_data_file = "./problem/GTOC12_Asteroids_Data.txt"
coeff_file = "./problem/bonus_coefficients_202306019.txt"

# Asteroid data
asteroid_df = DataFrame(CSV.File(asteroid_data_file; header=1, delim=' ', ignorerepeated=true))

# Add bonus coefficients and collected mass
coeff_df = DataFrame(CSV.File(coeff_file; header=false, delim=' ', ignorerepeated=true))
rename!(coeff_df, ["Bonus Coefficient", "Collected Mass (kg)"])
insertcols!(coeff_df, 1, :ID=>asteroid_df.ID) # Add IDs to data frame

# Join the two data frames to have a party
asteroid_df = ainnerjoin(asteroid_df, coeff_df, on=:ID)