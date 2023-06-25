using SPICE, AstroTime

export e⃗, ν, a, i, Ω, ω, RV2COE, COE2RV, body_osc_elt, M2EH, EH2ν, hyp_anom, ecc_anom, mean_anom #, get_planet_orbit

# Get osculating elements for a planet
function body_osc_elt(; planet::String, epoch::Epoch, frame=default_ref_frame, CB=default_CB_str)
    ET = Epoch_to_SPICE_ET(epoch)
    body_state, _ = spkgeo(planet, ET, frame, CB)
    return oscltx(body_state, ET, bodvrd(CB, "GM")[1])
end

function M2EH(; M, ecc, tol=eps(Float64))
    if ecc < 1
        EH = M + sign(π - mod(M, 2π)) * ecc
        EH₊ = EH # allocate outside loop scope
        while true # Vallado 4e Algorithm 2 p65
            EH₊ = EH + (M - EH + ecc * sin(EH)) / (1 - ecc * cos(EH))
            if tol > abs(EH₊ - EH)
                return EH₊
            end
            EH = EH₊
        end
    else
        EH = if ecc < 1.6
            M + sign(π - mod(M, 2π)) * ecc
        elseif (ecc < 3.6) & (abs(M) > π)
            M - sign(M) * ecc
        else
            M / (ecc - 1)
        end
        EH₊ = EH # allocate outside loop scope
        while true # Vallado 4e Algorithm 4 p71
            EH₊ = EH + (M - ecc * sinh(EH) + EH) / (ecc * cosh(EH) - 1)
            if tol > abs(EH₊ - EH)
                return EH₊
            end
            EH = EH₊
        end
    end
end

function EH2ν(; E_or_H, ecc)
    if ecc > 1
        return 2 * atan(sqrt((ecc + 1) / (ecc - 1)) * tanh(E_or_H / 2)) # Vallado 4e Eq. 2-36 (p56)
    else
        return 2 * atan(sqrt((1 + ecc) / (1 - ecc)) * tan(E_or_H / 2)) # Vallado 4e Eq. 2-13 (p48)
    end
end

function e⃗(; x⃗, μ_CB_or_CB_name)
    r⃗ = view(x⃗, 1:3)
    v⃗ = view(x⃗, 4:6)
    μ_CB = get_GM(μ_CB_or_CB_name)
    return ((norm(v⃗)^2 - μ_CB / norm(r⃗)) * r⃗ - (r⃗ ⋅ v⃗) * v⃗) / μ_CB # Vallado 4e Eq. 2-78 (p98)
end

function ν(; x⃗, μ_CB_or_CB_name)
    r⃗ = view(x⃗, 1:3)
    v⃗ = view(x⃗, 4:6)
    ecc = e⃗(x⃗=x⃗, μ_CB_or_CB_name=μ_CB_or_CB_name)
    ν̃ = acos((ecc ⋅ r⃗) / (norm(ecc) * norm(r⃗))) # Vallado 4e Eq. 2-86 (p100)
    return r⃗ ⋅ v⃗ > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace
end

function a(; x⃗, μ_CB_or_CB_name)
    r⃗ = view(x⃗, 1:3)
    v⃗ = view(x⃗, 4:6)
    return 1 / (2 / norm(r⃗) - norm(v⃗)^2 / get_GM(μ_CB_or_CB_name)) # Vallado 4e Eq. 2-74 (p96)
end

function i(; x⃗)
    r⃗ = view(x⃗, 1:3)
    v⃗ = view(x⃗, 4:6)
    h⃗ = r⃗ × v⃗
    return acos(normalize(h⃗)[3]) # Vallado 4e Eq. 2-82 (p99)
end

function Ω(; x⃗)
    r⃗ = view(x⃗, 1:3)
    v⃗ = view(x⃗, 4:6)
    h⃗ = r⃗ × v⃗
    n⃗ = [0, 0, 1] × h⃗ # Vallado 4e Eq. 2-83 (p99)
    RAAN = acos(normalize(n⃗)[1]) # Vallado 4e Eq. 2-84 (p99)
    return n⃗[2] > 0 ? RAAN : 2π - RAAN
end

function ω(; x⃗, μ_CB_or_CB_name)
    r⃗ = view(x⃗, 1:3)
    v⃗ = view(x⃗, 4:6)
    h⃗ = r⃗ × v⃗
    n⃗ = [0, 0, 1] × h⃗ # Vallado 4e Eq. 2-83 (p99)
    ecc = e⃗(r⃗=r⃗, v⃗=v⃗, μ_CB_or_CB_name=μ_CB_or_CB_name)
    AOP = acos((n⃗ ⋅ ecc) / (norm(n⃗) * norm(ecc))) # Vallado 4e Eq. 2-85 (p100)
    return ecc[3] > 0 ? AOP : 2π - AOP
end

"   COE = [a,e,i,Ω,ω,ν]"
function RV2COE(; x⃗, μ_CB_or_CB_name) # Vallado 4e Algorithm 9 (p113)
    r⃗ = view(x⃗, 1:3)
    v⃗ = view(x⃗, 4:6)
    r = norm(r⃗)
    v = norm(v⃗)
    h⃗ = r⃗ × v⃗
    n⃗ = [0, 0, 1] × h⃗ # Vallado 4e Eq. 2-83 (p99)
    μ_CB = get_GM(μ_CB_or_CB_name)
    sma = 1 / (2 / r - v^2 / μ_CB) # Vallado 4e Eq. 2-74 (p96)
    ecc = ((v^2 - μ_CB / r) * r⃗ - (r⃗ ⋅ v⃗) * v⃗) / μ_CB # Vallado 4e Eq. 2-78 (p98)
    e = norm(ecc)
    inc = acos(normalize(h⃗)[3]) # Vallado 4e Eq. 2-82 (p99)
    RAAN = acos(normalize(n⃗)[1]) # Vallado 4e Eq. 2-84 (p99)
    RAAN2 = n⃗[2] > 0 ? RAAN : 2π - RAAN
    AOP = acos((n⃗ ⋅ ecc) / (norm(n⃗) * e)) # Vallado 4e Eq. 2-85 (p100)
    AOP2 = ecc[3] > 0 ? AOP : 2π - AOP
    trueanom = acos((ecc ⋅ r⃗) / (e * r)) # Vallado 4e Eq. 2-86 (p100)
    trueanom2 = r⃗ ⋅ v⃗ > 0 ? trueanom : 2π - trueanom
    return [sma, e, inc, RAAN2, AOP2, trueanom2]
end

"   COE = [a,e,i,Ω,ω,ν]"
function COE2RV(; COE, μ_CB_or_CB_name) # Vallado 4e Algorithm 10 (p118)
    a, e, i, Ω, ω, ν = COE
    μ_CB = get_GM(μ_CB_or_CB_name)
    p = a * (1 - e^2)
    ci = cos(i)
    si = sin(i)
    cΩ = cos(Ω)
    sΩ = sin(Ω)
    cω = cos(ω)
    sω = sin(ω)
    cν = cos(ν)
    sν = sin(ν)
    f1 = (1 + e * cν)
    r̃ = [p * cν / f1, p * sν / f1, 0]
    f2 = sqrt(μ_CB / p)
    ṽ = [-f2 * sν, f2 * (e + cν), 0]
    IJK_PQW = [
        cΩ*cω-sΩ*sω*ci -cΩ*sω-sΩ*cω*ci sΩ*si
        sΩ*cω+cΩ*sω*ci -sΩ*sω+cΩ*cω*ci -cΩ*si
        sω*si cω*si ci
    ]
    x⃗ = Vector{Float64}(undef, 6)
    x⃗[1:3] = IJK_PQW * r̃
    x⃗[4:6] = IJK_PQW * ṽ
    return x⃗
end

function hyp_anom(; x⃗, μ_CB_or_CB_name)
    r⃗ = view(x⃗, 1:3)
    v⃗ = view(x⃗, 4:6)
    ecc = e⃗(x⃗=x⃗, μ_CB_or_CB_name=μ_CB_or_CB_name)
    e = norm(ecc)
    ν̃ = acos((ecc ⋅ r⃗) / (e * norm(r⃗))) # Vallado 4e Eq. 2-86 (p100)
    trueanom = r⃗ ⋅ v⃗ > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace
    return 2 * atanh(sqrt((e - 1) / (e + 1)) * tan(trueanom / 2)) # Vallado 4e Eq. 2-35 (p56)
end

function ecc_anom(; x⃗, μ_CB_or_CB_name)
    r⃗ = view(x⃗, 1:3)
    v⃗ = view(x⃗, 4:6)
    ecc = e⃗(x⃗=x⃗, μ_CB_or_CB_name=μ_CB_or_CB_name)
    e = norm(ecc)
    ν̃ = acos((ecc ⋅ r⃗) / (e * norm(r⃗))) # Vallado 4e Eq. 2-86 (p100)
    trueanom = r⃗ ⋅ v⃗ > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace
    return 2 * atan(sqrt((1 - e) / (1 + e)) * tan(trueanom / 2)) # Vallado 4e Eq. 2-14 (p48)
end

function mean_anom(; x⃗, μ_CB_or_CB_name)
    r⃗ = view(x⃗, 1:3)
    v⃗ = view(x⃗, 4:6)
    ecc = e⃗(x⃗=x⃗, μ_CB_or_CB_name=μ_CB_or_CB_name)
    e = norm(ecc)
    ν̃ = acos((ecc ⋅ r⃗) / (e * norm(r⃗))) # Vallado 4e Eq. 2-86 (p100)
    trueanom = r⃗ ⋅ v⃗ > 0 ? ν̃ : 2π - ν̃ # Correct for halfspace
    if e < 1
        E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(trueanom / 2)) # Vallado 4e Eq. 2-14 (p48)
        return E - e * sin(E) # Vallado 4e Eq. 2-4 (p45)
    elseif e > 1
        H = 2 * atanh(sqrt((e - 1) / (e + 1)) * tan(trueanom / 2)) # Vallado 4e Eq. 2-35 (p56)
        return e * sinh(H) - H # Vallado 4e Eq. 2-38 (p57)
    end
end