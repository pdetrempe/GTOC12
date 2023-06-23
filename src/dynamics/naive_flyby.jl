using LinearAlgebra, SPICE

export e⃗, ν, a, i, Ω, ω, RV2COE, COE2RV, body_SOI, hyp_anom, ecc_anom, mean_anom, hyp_turn_angle, hyp_periapsis, hyp_exit_r⃗, hyp_exit_v⃗, hyp_exit_x⃗, flyby_TOF, naive_flyby

get_GM(μ_CB_or_CB_name) = typeof(μ_CB_or_CB_name) != String ? μ_CB_or_CB_name : bodvrd(μ_CB_or_CB_name,"GM")[1] # if GM provided directly (is a number), use it, else retrieve from body name (String)

function e⃗(;x⃗,μ_CB_or_CB_name)
    r⃗ = view(x⃗,1:3)
    v⃗ = view(x⃗,4:6)
    μ_CB = get_GM(μ_CB_or_CB_name)
    return ((norm(v⃗)^2 - μ_CB/norm(r⃗))*r⃗ - (r⃗⋅v⃗)*v⃗)/μ_CB # Vallado 4e Eq. 2-78 (p98)
end

function ν(;x⃗,μ_CB_or_CB_name)
    r⃗ = view(x⃗,1:3)
    v⃗ = view(x⃗,4:6)
    ecc = e⃗(x⃗=x⃗,μ_CB_or_CB_name=μ_CB_or_CB_name)
    ν̃ = acos((ecc⋅r⃗)/(norm(ecc)*norm(r⃗))) # Vallado 4e Eq. 2-86 (p100)
    return r⃗⋅v⃗ > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace
end

function a(;x⃗,μ_CB_or_CB_name)
    r⃗ = view(x⃗,1:3)
    v⃗ = view(x⃗,4:6)
    return 1/(2/norm(r⃗) - norm(v⃗)^2/get_GM(μ_CB_or_CB_name)) # Vallado 4e Eq. 2-74 (p96)
end

function i(;x⃗)
    r⃗ = view(x⃗,1:3)
    v⃗ = view(x⃗,4:6)
    h⃗ = r⃗ × v⃗
    return acos(normalize(h⃗)[3]) # Vallado 4e Eq. 2-82 (p99)
end

function Ω(;x⃗)
    r⃗ = view(x⃗,1:3)
    v⃗ = view(x⃗,4:6)
    h⃗ = r⃗ × v⃗
    n⃗ = [0,0,1] × h⃗ # Vallado 4e Eq. 2-83 (p99)
    RAAN = acos(normalize(n⃗)[1]) # Vallado 4e Eq. 2-84 (p99)
    return n⃗[2] > 0 ? RAAN : 2π - RAAN
end

function ω(;x⃗,μ_CB_or_CB_name)
    r⃗ = view(x⃗,1:3)
    v⃗ = view(x⃗,4:6)
    h⃗ = r⃗ × v⃗
    n⃗ = [0,0,1] × h⃗ # Vallado 4e Eq. 2-83 (p99)
    ecc = e⃗(r⃗=r⃗,v⃗=v⃗,μ_CB_or_CB_name=μ_CB_or_CB_name)
    AOP = acos((n⃗⋅ecc)/(norm(n⃗)*norm(ecc))) # Vallado 4e Eq. 2-85 (p100)
    return ecc[3] > 0 ? AOP : 2π - AOP
end

"   COE = [a,e,i,Ω,ω,ν]"
function RV2COE(;x⃗,μ_CB_or_CB_name) # Vallado 4e Algorithm 9 (p113)
    r⃗ = view(x⃗,1:3)
    v⃗ = view(x⃗,4:6)
    r = norm(r⃗)
    v = norm(v⃗)
    h⃗ = r⃗ × v⃗
    n⃗ = [0,0,1] × h⃗ # Vallado 4e Eq. 2-83 (p99)
    μ_CB = get_GM(μ_CB_or_CB_name)
    sma = 1/(2/r - v^2/μ_CB) # Vallado 4e Eq. 2-74 (p96)
    ecc = ((v^2 - μ_CB/r)*r⃗ - (r⃗⋅v⃗)*v⃗)/μ_CB # Vallado 4e Eq. 2-78 (p98)
    e = norm(ecc)
    inc = acos(normalize(h⃗)[3]) # Vallado 4e Eq. 2-82 (p99)
    RAAN = acos(normalize(n⃗)[1]) # Vallado 4e Eq. 2-84 (p99)
    RAAN2 = n⃗[2] > 0 ? RAAN : 2π - RAAN
    AOP = acos((n⃗⋅ecc)/(norm(n⃗)*e)) # Vallado 4e Eq. 2-85 (p100)
    AOP2 = ecc[3] > 0 ? AOP : 2π - AOP
    trueanom = acos((ecc⋅r⃗)/(e*r)) # Vallado 4e Eq. 2-86 (p100)
    trueanom2 = r⃗⋅v⃗ > 0 ? trueanom : 2π - trueanom
    return [sma,e,inc,RAAN2,AOP2,trueanom2]
end

"   COE = [a,e,i,Ω,ω,ν]"
function COE2RV(;COE,μ_CB_or_CB_name) # Vallado 4e Algorithm 10 (p118)
    a,e,i,Ω,ω,ν = COE
    μ_CB = get_GM(μ_CB_or_CB_name)
    p = a*(1-e^2)
    ci = cos(i); si = sin(i)
    cΩ = cos(Ω); sΩ = sin(Ω)
    cω = cos(ω); sω = sin(ω)
    cν = cos(ν); sν = sin(ν)
    f1 = (1+e*cν)
    r̃ = [p*cν/f1, p*sν/f1, 0]
    f2 = sqrt(μ_CB/p)
    ṽ = [-f2*sν, f2*(e+cν), 0]
    IJK_PQW = [
        cΩ*cω-sΩ*sω*ci  -cΩ*sω-sΩ*cω*ci sΩ*si
        sΩ*cω+cΩ*sω*ci  -sΩ*sω+cΩ*cω*ci -cΩ*si
        sω*si           cω*si           ci
    ]
    x⃗ = Vector{Float64}(undef,6)
    x⃗[1:3] = IJK_PQW * r̃
    x⃗[4:6] = IJK_PQW * ṽ
    return x⃗
end

function body_SOI(;CB::String,orbiting_body::String)
    orbiting_body_state = spkgeo(bodn2c(orbiting_body),0,base_ref_frame,bodn2c(CB))[1]
    orbiting_body_GM = bodvrd(orbiting_body,"GM")[1]
    CB_GM = bodvrd(CB,"GM")[1]
    orbiting_body_a = 1/(2/norm(orbiting_body_state[1:3]) - norm(orbiting_body_state[4:6])^2/CB_GM) # Vallado 4e Eq. 2-74 (p96)
    return orbiting_body_a*(orbiting_body_GM/CB_GM)^0.4 # Vallado 4e Eq. 12-2 (p948)
end

function hyp_anom(;x⃗,μ_CB_or_CB_name)
    r⃗ = view(x⃗,1:3)
    v⃗ = view(x⃗,4:6)
    ecc = e⃗(x⃗=x⃗,μ_CB_or_CB_name=μ_CB_or_CB_name)
    e = norm(ecc)
    ν̃ = acos((ecc⋅r⃗)/(e*norm(r⃗))) # Vallado 4e Eq. 2-86 (p100)
    trueanom = r⃗⋅v⃗ > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace
    return 2*atanh(sqrt((e-1)/(e+1)) * tan(trueanom/2)) # Vallado 4e Eq. 2-35 (p56)
end

function ecc_anom(;x⃗,μ_CB_or_CB_name)
    r⃗ = view(x⃗,1:3)
    v⃗ = view(x⃗,4:6)
    ecc = e⃗(x⃗=x⃗,μ_CB_or_CB_name=μ_CB_or_CB_name)
    e = norm(ecc)
    ν̃ = acos((ecc⋅r⃗)/(e*norm(r⃗))) # Vallado 4e Eq. 2-86 (p100)
    trueanom = r⃗⋅v⃗ > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace
    return 2*atan(sqrt((1-e)/(1+e)) * tan(trueanom/2)) # Vallado 4e Eq. 2-14 (p48)
end

function mean_anom(;x⃗,μ_CB_or_CB_name)
    r⃗ = view(x⃗,1:3)
    v⃗ = view(x⃗,4:6)
    ecc = e⃗(x⃗=x⃗,μ_CB_or_CB_name=μ_CB_or_CB_name)
    e = norm(ecc)
    ν̃ = acos((ecc⋅r⃗)/(e*norm(r⃗))) # Vallado 4e Eq. 2-86 (p100)
    trueanom = r⃗⋅v⃗ > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace
    if e < 1 
        E = 2*atan(sqrt((1-e)/(1+e)) * tan(trueanom/2)) # Vallado 4e Eq. 2-14 (p48)
        return E - e*sin(E) # Vallado 4e Eq. 2-4 (p45)
    elseif e > 1
        H = 2*atanh(sqrt((e-1)/(e+1)) * tan(trueanom/2)) # Vallado 4e Eq. 2-35 (p56)
        return e*sinh(H) - H # Vallado 4e Eq. 2-38 (p57)
    end
end

function hyp_turn_angle(;x⃗,μ_CB_or_CB_name)
    ecc = e⃗(x⃗=x⃗,μ_CB_or_CB_name=μ_CB_or_CB_name)
    return 2*asin(1/norm(ecc)) # Vallado 4e Eq. 2-28 (p53)
end

function hyp_periapsis(;x⃗∞,μ_CB_or_CB_name)
    v⃗∞ = view(x⃗∞,4:6)
    ecc = e⃗(x⃗=x⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name)
    turn_angle = 2*asin(1/norm(ecc)) # Vallado 4e Eq. 2-28 (p53)
    return get_GM(μ_CB_or_CB_name)/norm(v⃗∞)^2 * (1/cos((π-turn_angle)/2)-1) # Vallado 4e Eq. 12-12 (p959)
end

# h⃗(;r⃗,v⃗) = r⃗ × v⃗ # Specific angular momentum
# ĥ(;r⃗,v⃗) = normalize(h⃗(r⃗=r⃗,v⃗=v⃗)) # Specific angular momentum unit vector

function hyp_exit_r⃗(;x⃗∞,μ_CB_or_CB_name)
    r⃗∞ = view(x⃗∞,1:3)
    v⃗∞ = view(x⃗∞,4:6)
    return vrotv( # Mirror the radius vector "at infinity" about the periapsis vector
        r⃗∞, # To be rotated
        e⃗(x⃗=x⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name), # periapsis vector (not required to be unit vector)
        π # 180 degrees
    )
end

function hyp_exit_v⃗(;x⃗∞,μ_CB_or_CB_name)
    r⃗∞ = view(x⃗∞,1:3)
    v⃗∞ = view(x⃗∞,4:6)
    return vrotv(
        v⃗∞,
        r⃗∞ × v⃗∞, # axis of rotation as specific angular momentum vector
        hyp_turn_angle(x⃗=x⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name) # turn by the hyperbolic turn angle
    )
end

function hyp_exit_x⃗(;x⃗∞,μ_CB_or_CB_name)
    r⃗∞ = view(x⃗∞,1:3)
    v⃗∞ = view(x⃗∞,4:6)
    ecc = e⃗(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name)
    exit_x⃗ = Vector{Float64}(undef,6)
    exit_x⃗[1:3] = vrotv( # Mirror the radius vector "at infinity" about the periapsis vector
    r⃗∞, # To be rotated
    ecc, # periapsis vector (not required to be unit vector)
    π # 180 degrees
    )
    exit_x⃗[4:6] = vrotv(
        v⃗∞,
        r⃗∞ × v⃗∞, # axis of rotation as specific angular momentum vector
        hyp_turn_angle(e=norm(ecc)) # turn by the hyperbolic turn angle
    )
    return exit_x⃗
end

function flyby_TOF(;x⃗∞,μ_CB_or_CB_name)
    r⃗∞ = view(x⃗∞,1:3)
    v⃗∞ = view(x⃗∞,4:6)
    μ_CB = get_GM(μ_CB_or_CB_name)
    ecc = ((norm(v⃗∞)^2 - μ_CB/norm(r⃗∞))*r⃗ - (r⃗∞⋅v⃗∞)*v⃗∞)/μ_CB # Vallado 4e Eq. 2-78 (p98)
    e = norm(ecc)
    ν̃ = acos((ecc⋅r⃗∞)/(e*norm(r⃗∞))) # Vallado 4e Eq. 2-86 (p100)
    νin = r⃗∞⋅v⃗∞ > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace; true anomaly at SOI entry
    Hin = 2*atanh(sqrt((e-1)/(e+1)) * tan(νin/2)) # Hyperbolic anomaly at SOI entry
    sma = 1/(2/norm(r⃗∞) - norm(v⃗∞)^2/μ_CB) # Vallado 4e Eq. 2-74 (p96)
    return sqrt(-sma^3/μ_CB)*2*(Hin - e*sinh(Hin)) # Vallado 4e Eq. 2-39 (p57) with simplifications since Hout = -Hin
end

function naive_flyby(;x⃗_inrt,epoch_et,flyby_body,CB=default_CB_str,inrt_frame=default_ref_frame)
    fbbdyc = bodn2c(flyby_body)
    CBc = bodn2c(CB)

    x⃗_fbbdy = spkgeo(fbbdyc,epoch_et,inrt_frame,CBc)[1]
    x⃗∞in = x⃗_inrt - x⃗_fbbdy
    r⃗∞ = view(x⃗∞in,1:3)
    r∞ = norm(r⃗∞)
    v⃗∞ = view(x⃗∞in,4:6)
    v∞2 = norm(v⃗∞)^2
    rdv = r⃗∞⋅v⃗∞
    μ_fbbdy = bodvrd(flyby_body,"GM")[1]
    ecc = ((v∞2 - μ_fbbdy/r∞)*r⃗∞ - (rdv)*v⃗∞)/μ_fbbdy # Vallado 4e Eq. 2-78 (p98)
    e = norm(ecc)
    
    ν̃ = acos((ecc⋅r⃗∞)/(e*r∞)) # Vallado 4e Eq. 2-86 (p100)
    νin = rdv > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace; true anomaly at SOI entry
    Hin = 2*atanh(sqrt((e-1)/(e+1)) * tan(νin/2)) # Hyperbolic anomaly at SOI entry
    sma = 1/(2/r∞ - v∞2/μ_fbbdy) # Vallado 4e Eq. 2-74 (p96)
    Δt = sqrt(-sma^3/μ_CB)*2*(Hin - e*sinh(Hin)) # Vallado 4e Eq. 2-39 (p57) with simplifications since Hout = -Hin
    epoch_et_out = epoch_et + Δt

    x⃗∞out = Vector{Float64}(undef,6)
    x⃗∞out[1:3] = vrotv( # Mirror the radius vector "at infinity" about the periapsis vector
    r⃗∞, # To be rotated
    ecc, # periapsis vector (not required to be unit vector)
    π # 180 degrees
    )
    x⃗∞out[4:6] = vrotv(
        v⃗∞,
        r⃗∞ × v⃗∞, # axis of rotation as specific angular momentum vector
        2*asin(1/e) # Vallado 4e Eq. 2-28 (p53), turn by the hyperbolic turn angle
    )
    x⃗_fbbdy_out = spkgeo(fbbdyc,epoch_et_out,inrt_frame,CBc)[1]

    return x⃗∞out + x⃗_fbbdy_out, epoch_et_out
end