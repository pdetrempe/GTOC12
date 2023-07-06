using SPICE

function get_planet_df()
    planet_df = DataFrame([
        Pair(itr...)
        for itr in zip(
        [:ID, :NAME, :sma, :ecc, :inc, :LAN, :argperi, :true_anom, :mean_anom],
        [Int[],String[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[]]
    )]
    )
    plts = [1,2,3,4,5,6,7,8,9,199,299,399]
    for idx in eachindex(plts)
        x_plt=get_planet_state(plts[idx],ET₀)
        plt_oe = RV2COE(x⃗=x_plt)
        plt_M = mean_anom(x⃗=x_plt)
        push!(planet_df,[plts[idx], bodc2s(plts[idx]), plt_oe..., plt_M])
    end

    planet_df[!, :MJD] .= 64328
    planet_df[!, :TT] .= TT_start
    planet_df[!, :ET] .= str2et(TT_start)

    return planet_df
end