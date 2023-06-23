using LinearAlgebra, UnPack

export keplerian2MEE, MME2keplerian, MEE2Cartesian, Cartesian2MEE, EOM_MEE!

## Intermediate quantities used in MEE calculations
get_q(; f, g, L) = 1 + f * cos(L) + g * sin(L)
get_s(; h, k) = sqrt(1 + h^2 + k^2)
get_α(; h, k) = h^2 - k^2
get_r(; p, w) = p / w
get_w(; f, g, L) = 1 + f * cos(L) + g * sin(L)

## Orbital element conversions
function keplerian2MEE(; a, e, i, Ω, ω, ν)
    """
    Converts from classical Keplerian elements to modified equinoctial elements
    """
    p = a * (1 - e^2)
    f = e * cos(Ω + ω)
    g = e * sin(Ω + ω)
    h = tan(i / 2) * cos(Ω)
    k = tan(i / 2) * sin(Ω)
    l = Ω + ω + ν

    [p, f, g, h, k, l]

end

function MEE2keplerian(; p, f, g, h, k, l)
    """
    Converts from modified equinoctial elements to classical Keplerian 
    See https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
    """
    a = p / (1 - f^2 - g^2)
    e = sqrt(f^2 + g^2)
    i = atan(2 * sqrt(h^2 + k^2), 1 - h^2 - k^2)
    Ω = atan(k, h)
    ω = atan(g * h - f * k, f * h + g * k)
    ν = l - (Ω + ω)



    [a, e, i, Ω, ω, ν]
end

function MEE2Cartesian(MEE; μ)
    """
    Converts from modified equinoctial elements to Cartesian coordinates
    See https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
    Equations 3a & 3b
    https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit/blob/master/src/modified_equinoctial_module.f90
    """    
    p, f, g, h, k, L = MEE

    fhat= zeros(eltype(k), 3)
    ghat = zeros(eltype(k), 3)

    kk      = k^2
    hh      = h^2
    tkh     = 2*k*h
    s2      = 1 + hh + kk
    cL      = cos(L)
    sL      = sin(L)
    w       = 1 + f*cL + g*sL
    r       = p/w
    smp     = sqrt(μ/p)
    fhat[1] = 1-kk+hh
    fhat[3] = -2*k
    ghat[1] = tkh
    ghat[2] = 1+kk-hh
    fhat[2] = tkh
    ghat[3] = 2*h
    fhat    = fhat/s2
    ghat    = ghat/s2
    x       = r*cL
    y       = r*sL
    xdot    = -smp * (g + sL)
    ydot    =  smp * (f + cL)

    r⃗ = x*fhat + y*ghat
    v⃗ = xdot*fhat + ydot*ghat
    
    vcat(r⃗, v⃗)
end

function Cartesian2MEE(x⃗; μ)
    """
    https://degenerateconic.com/modified-equinoctial-elements.html
    """
    r = x⃗[1:3]
    v = x⃗[4:6]

    fhat= zeros(eltype(r), 3)
    ghat = zeros(eltype(r), 3)

    rdv      = dot(r,v)
    rmag     = norm(r)
    rhat     = r/rmag
    hvec     = cross(r,v)
    hmag     = norm(hvec)
    hhat     = hvec/hmag
    vhat     = (rmag*v - rdv*rhat)/hmag
    p        = hmag*hmag / μ
    k        = hhat[1]/(1 + hhat[3])
    h        = -hhat[2]/(1 + hhat[3])
    kk       = k*k
    hh       = h*h
    s2       = 1+hh+kk
    tkh      = 2*k*h
    ecc      = cross(v,hvec)/μ - rhat
    fhat[1]  = 1-kk+hh
    fhat[2]  = tkh
    fhat[3]  = -2*k
    ghat[1]  = tkh
    ghat[2]  = 1+kk-hh
    ghat[3]  = 2*h
    fhat     = fhat/s2
    ghat     = ghat/s2
    f        = dot(ecc,fhat)
    g        = dot(ecc,ghat)
    L        = mod2pi(atan(rhat[2]-vhat[1],rhat[1]+vhat[2]))

    return [p;f;g;h;k;L]

end


## MEE dynamics

function A_equinoctial(MEE; μ)
    """
    See: https://ai.jpl.nasa.gov/public/documents/papers/AAS-22-015-Paper.pdf
    Eq 2
    """

    p, f, g, h, k, l = MEE
    q = get_q(f=f, g=g, L=l)

    A = [0; 0; 0; 0; 0; sqrt(μ * p) * (q / p)^2]
end

function B_equinoctial(MEE; μ)
    """
    See: https://ai.jpl.nasa.gov/public/documents/papers/AAS-22-015-Paper.pdf
    Eq 3
    """

    p, f, g, h, k, l = MEE
    q = get_q(f=f, g=g, L=l)
    s = get_s(h=h, k=k)

    B = [0 2*p/q*sqrt(p / μ) 0
        sqrt(p / μ)*sin(l) sqrt(p / μ)*q*((q+1)*cos(l)+f) -sqrt(p / μ)*g/q*(h*sin(l)-k*cos(l))
        -sqrt(p / μ)*cos(l) sqrt(p / μ)*q*((q+1)*sin(l)+g) sqrt(p / μ)*f/q*(h*sin(l)-k*cos(l))
        0 0 sqrt(p / μ)*s*cos(l)/(2*q)
        0 0 sqrt(p / μ)*s*sin(l)/(2*q)
        0 0 sqrt(p / μ)*1/q*(h*sin(l)-k*cos(l))]
end


function tangential_firing(MEE; params) # Get control thrust direction and magnitude [0, 1]
    @unpack μ = params # problem parameters

    x⃗ = MEE2Cartesian(MEE; μ)
    v⃗ = view(x⃗, 4:6)
    R_inrt2lvlh = DCM_inertial_to_lvlh(x⃗)

    # Just roll with tangential firing to sanity check
    û_LVLH = R_inrt2lvlh * v⃗ / norm(v⃗) # thrust along velocity vector (in LVLH frame)
    δ = 1         # Full throttle

    û_LVLH, δ
end

function EOM_MEE!(ẋ, x, p, t)
    """
    See: https://ai.jpl.nasa.gov/public/documents/papers/AAS-22-015-Paper.pdf
    Eq 1
    """
    @unpack μ, c, T_max = p # problem parameters
    MEE = x[1:6]
    m = x[7]

    A = A_equinoctial(MEE; μ=μ)
    B = B_equinoctial(MEE; μ=μ)
    û, δ = get_control(MEE; params=p) # Get control thrust direction and magnitude [0, 1]

    MEE_dot = A + T * δ / m * B * û
    ṁ = -δ * T / (c)

    ẋ[:] = vcat(MEE_dot, ṁ)

end