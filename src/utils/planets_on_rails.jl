using SPICE

function get_planet_df(plts = [1,2,3,4,5,6,7,8,9,199,299,399])
    planet_df = DataFrame([
        Pair(itr...)
        for itr in zip(
        [:ID, :NAME, :sma, :ecc, :inc, :LAN, :argperi, :true_anom, :mean_anom, :SOI],
        [Int[],String[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[]]
    )]
    )
    for idx in eachindex(plts)
        pc = plts[idx]
        pn = bodc2s(plts[idx])
        x_plt=get_planet_state(pc,ET₀)
        plt_oe = RV2COE(x⃗=x_plt)
        plt_M = mean_anom(x⃗=x_plt)
        plt_SOI = body_SOI(CB=default_CB_str,orbiting_body=pn)
        push!(planet_df,[pc, pn, plt_oe..., plt_M, plt_SOI])
    end

    planet_df[!, :MJD] .= 64328
    planet_df[!, :TT] .= TT_start
    planet_df[!, :ET] .= str2et(TT_start)

    return planet_df
end