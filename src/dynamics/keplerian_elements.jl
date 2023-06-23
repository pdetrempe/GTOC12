using SPICE, AstroTime

export body_osc_elt, M2EH, EH2ν #, get_planet_orbit

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
