
export body_SOI, hyp_turn_angle, hyp_periapsis, hyp_exit_r⃗, hyp_exit_v⃗, hyp_exit_x⃗, flyby_TOF, naive_flyby

get_GM(μ_CB_or_CB_name) = typeof(μ_CB_or_CB_name) != String ? μ_CB_or_CB_name : bodvrd(μ_CB_or_CB_name, "GM")[1] # if GM provided directly (is a number), use it, else retrieve from body name (String)

function body_SOI(; CB::String, orbiting_body::String)
    orbiting_body_state = spkgeo(bodn2c(orbiting_body), 0, default_ref_frame, bodn2c(CB))[1]
    orbiting_body_GM = bodvrd(orbiting_body, "GM")[1]
    CB_GM = bodvrd(CB, "GM")[1]
    orbiting_body_a = 1 / (2 / norm(orbiting_body_state[1:3]) - norm(orbiting_body_state[4:6])^2 / CB_GM) # Vallado 4e Eq. 2-74 (p96)
    return orbiting_body_a * (orbiting_body_GM / CB_GM)^0.4 # Vallado 4e Eq. 12-2 (p948)
end

function hyp_turn_angle(; x⃗, μ_CB_or_CB_name)
    ecc = e⃗(x⃗=x⃗, μ_CB_or_CB_name=μ_CB_or_CB_name)
    return 2 * asin(1 / norm(ecc)) # Vallado 4e Eq. 2-28 (p53)
end

function hyp_periapsis(; x⃗∞, μ_CB_or_CB_name)
    v⃗∞ = view(x⃗∞, 4:6)
    ecc = e⃗(x⃗=x⃗∞, μ_CB_or_CB_name=μ_CB_or_CB_name)
    turn_angle = 2 * asin(1 / norm(ecc)) # Vallado 4e Eq. 2-28 (p53)
    return get_GM(μ_CB_or_CB_name) / norm(v⃗∞)^2 * (1 / cos((π - turn_angle) / 2) - 1) # Vallado 4e Eq. 12-12 (p959)
end

# h⃗(;r⃗,v⃗) = r⃗ × v⃗ # Specific angular momentum
# ĥ(;r⃗,v⃗) = normalize(h⃗(r⃗=r⃗,v⃗=v⃗)) # Specific angular momentum unit vector

function hyp_exit_r⃗(; x⃗∞, μ_CB_or_CB_name)
    r⃗∞ = view(x⃗∞, 1:3)
    return vrotv( # Mirror the radius vector "at infinity" about the periapsis vector
        r⃗∞, # To be rotated
        e⃗(x⃗=x⃗∞, μ_CB_or_CB_name=μ_CB_or_CB_name), # periapsis vector (not required to be unit vector)
        π # 180 degrees
    )
end

function hyp_exit_v⃗(; x⃗∞, μ_CB_or_CB_name)
    r⃗∞ = view(x⃗∞, 1:3)
    v⃗∞ = view(x⃗∞, 4:6)
    return vrotv(
        v⃗∞,
        r⃗∞ × v⃗∞, # axis of rotation as specific angular momentum vector
        hyp_turn_angle(x⃗=x⃗∞, μ_CB_or_CB_name=μ_CB_or_CB_name) # turn by the hyperbolic turn angle
    )
end

function hyp_exit_x⃗(; x⃗∞, μ_CB_or_CB_name)
    r⃗∞ = view(x⃗∞, 1:3)
    v⃗∞ = view(x⃗∞, 4:6)
    ecc = e⃗(x⃗=x⃗∞, μ_CB_or_CB_name=μ_CB_or_CB_name)
    exit_x⃗ = Vector{Float64}(undef, 6)
    exit_x⃗[1:3] = vrotv( # Mirror the radius vector "at infinity" about the periapsis vector
        r⃗∞, # To be rotated
        ecc, # periapsis vector (not required to be unit vector)
        π # 180 degrees
    )
    exit_x⃗[4:6] = vrotv(
        v⃗∞,
        r⃗∞ × v⃗∞, # axis of rotation as specific angular momentum vector
        2 * asin(1 / norm(ecc)) # Vallado 4e Eq. 2-28 (p53) # turn by the hyperbolic turn angle
    )
    return exit_x⃗
end

function flyby_TOF(; x⃗∞, μ_CB_or_CB_name)
    r⃗∞ = view(x⃗∞, 1:3)
    v⃗∞ = view(x⃗∞, 4:6)
    μ_CB = get_GM(μ_CB_or_CB_name)
    ecc = ((norm(v⃗∞)^2 - μ_CB / norm(r⃗∞)) * r⃗ - (r⃗∞ ⋅ v⃗∞) * v⃗∞) / μ_CB # Vallado 4e Eq. 2-78 (p98)
    e = norm(ecc)
    ν̃ = acos((ecc ⋅ r⃗∞) / (e * norm(r⃗∞))) # Vallado 4e Eq. 2-86 (p100)
    νin = r⃗∞ ⋅ v⃗∞ > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace; true anomaly at SOI entry
    Hin = 2 * atanh(sqrt((e - 1) / (e + 1)) * tan(νin / 2)) # Hyperbolic anomaly at SOI entry
    sma = 1 / (2 / norm(r⃗∞) - norm(v⃗∞)^2 / μ_CB) # Vallado 4e Eq. 2-74 (p96)
    return sqrt(-sma^3 / μ_CB) * 2 * (Hin - e * sinh(Hin)) # Vallado 4e Eq. 2-39 (p57) with simplifications since Hout = -Hin
end

function naive_flyby(; x⃗_inrt, epoch_et, flyby_body::String, CB::String=default_CB_str, inrt_frame=default_ref_frame)
    fbbdyc = bodn2c(flyby_body)
    CBc = bodn2c(CB)

    x⃗_fbbdy = spkgeo(fbbdyc, epoch_et, inrt_frame, CBc)[1]
    x⃗∞in = x⃗_inrt - x⃗_fbbdy
    r⃗∞ = view(x⃗∞in, 1:3)
    r∞ = norm(r⃗∞)
    v⃗∞ = view(x⃗∞in, 4:6)
    v∞2 = norm(v⃗∞)^2
    rdv = r⃗∞ ⋅ v⃗∞
    μ_fbbdy = bodvrd(flyby_body, "GM")[1]
    ecc = ((v∞2 - μ_fbbdy / r∞) * r⃗∞ - (rdv) * v⃗∞) / μ_fbbdy # Vallado 4e Eq. 2-78 (p98)
    e = norm(ecc)

    ν̃ = acos((ecc ⋅ r⃗∞) / (e * r∞)) # Vallado 4e Eq. 2-86 (p100)
    νin = rdv > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace; true anomaly at SOI entry
    Hin = 2 * atanh(sqrt((e - 1) / (e + 1)) * tan(νin / 2)) # Hyperbolic anomaly at SOI entry
    sma = 1 / (2 / r∞ - v∞2 / μ_fbbdy) # Vallado 4e Eq. 2-74 (p96)
    Δt = sqrt(-sma^3 / μ_CB) * 2 * (Hin - e * sinh(Hin)) # Vallado 4e Eq. 2-39 (p57) with simplifications since Hout = -Hin
    epoch_et_out = epoch_et + Δt

    x⃗∞out = Vector{Float64}(undef, 6)
    x⃗∞out[1:3] = vrotv( # Mirror the radius vector "at infinity" about the periapsis vector
        r⃗∞, # To be rotated
        ecc, # periapsis vector (not required to be unit vector)
        π # 180 degrees
    )
    x⃗∞out[4:6] = vrotv(
        v⃗∞,
        r⃗∞ × v⃗∞, # axis of rotation as specific angular momentum vector
        2 * asin(1 / e) # Vallado 4e Eq. 2-28 (p53), turn by the hyperbolic turn angle
    )
    x⃗_fbbdy_out = spkgeo(fbbdyc, epoch_et_out, inrt_frame, CBc)[1]

    return x⃗∞out + x⃗_fbbdy_out, epoch_et_out
end


function naive_flyby(; x⃗_inrt, epoch_et, flyby_body::Planet, rp=nothing)
    # TODO, take in a target v∞_out, calculate turn angle, ϕ, use that to calculate rp
    if rp === nothing
        rp = flyby_body.r_peri_min
    end

    x⃗_fbbdy = get_body_state(flyby_body; ET=epoch_et) #spkgeo(fbbdyc, epoch_et, inrt_frame, CBc)[1]
    x⃗∞in = x⃗_inrt - x⃗_fbbdy
    r⃗∞ = view(x⃗∞in, 1:3)
    r∞ = norm(r⃗∞)
    v⃗∞ = view(x⃗∞in, 4:6)
    v∞2 = norm(v⃗∞)^2
    rdv = r⃗∞ ⋅ v⃗∞
    μ_fbbdy = flyby_body.μ #bodvrd(flyby_body, "GM")[1]
    ecc = ((v∞2 - μ_fbbdy / r∞) * r⃗∞ - (rdv) * v⃗∞) / μ_fbbdy # Vallado 4e Eq. 2-78 (p98)
    e = norm(ecc)

    # Use rp to calculate eccentricity
    e = 1 + v∞2*rp/μ_fbbdy
    println("e $e")



    # ν̃ = acos((ecc ⋅ r⃗∞) / (e * r∞)) # Vallado 4e Eq. 2-86 (p100)
    # νin = rdv > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace; true anomaly at SOI entry
    # Hin = 2 * atanh(sqrt((e - 1) / (e + 1)) * tan(νin / 2)) # Hyperbolic anomaly at SOI entry
    # sma = 1 / (2 / r∞ - v∞2 / μ_fbbdy) # Vallado 4e Eq. 2-74 (p96)
    # Δt = sqrt(-sma^3 / μ_CB) * 2 * (Hin - e * sinh(Hin)) # Vallado 4e Eq. 2-39 (p57) with simplifications since Hout = -Hin
    # epoch_et_out = epoch_et # + Δt

    x⃗∞out = Vector{Float64}(undef, 6)
    x⃗∞out[1:3] = vrotv( # Mirror the radius vector "at infinity" about the periapsis vector
        r⃗∞, # To be rotated
        ecc, # periapsis vector (not required to be unit vector)
        π # 180 degrees
    )
    x⃗∞out[4:6] = vrotv(
        v⃗∞,
        r⃗∞ × v⃗∞, # axis of rotation as specific angular momentum vector
        2 * asin(1 / e) # Vallado 4e Eq. 2-28 (p53), turn by the hyperbolic turn angle
    )
    # x⃗_fbbdy_out = spkgeo(fbbdyc, epoch_et_out, inrt_frame, CBc)[1]

    return x⃗∞out + x⃗_fbbdy
end