using LinearAlgebra

export keplerian2MEE, MME2keplerian, MEE2Cartesian, Cartesian2MEE, EOM_MEE!

## Intermediate quantities used in MEE calculations
get_q(;f, g, L) = 1 + f*cos(L) + g*sin(L)
get_s(;h, k) = sqrt(1 + h^2 + k^2)
get_α(;h, k) = h^2 - k^2
get_r(;p, w) = p/w
get_w(;f, g, L) = 1+f*cos(L)+g*sin(L)

## Orbital element conversions
function keplerian2MEE(;a, e, i, Ω, ω, ν)
    """
    Converts from classical Keplerian elements to modified equinoctial elements
    """
    p = a*(1-e^2)
    f = e*cos(Ω + ω)
    g = e*sin(Ω + ω)
    h = tan(i/2)*cos(Ω)
    k = tan(i/2)*sin(Ω)
    l = Ω + ω + ν

    [p, f, g, h, k, l]

end

function MEE2keplerian(;p, f, g, h, k, l)
    """
    Converts from modified equinoctial elements to classical Keplerian 
    See https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
    """
    a = p/(1 - f^2 - g^2)
    e = sqrt(f^2 + g^2)
    i = atan( 2*sqrt(h^2 + k^2), 1-h^2-k^2)
    Ω = atan( k,h )
    ω = atan( g*h-f*k, f*h+g*k )
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
    p, f, g, h, k, l = MEE
    
    α = get_α(h=h, k=k)
    s = get_s(h=h, k=k)
    w = get_w(;f=f, g=g, L=l) 

    r = p/w
    r⃗ = r/(s^2) * [cos(l) + α^2 * cos(l) + 2*h*k*sin(l)
                    sin(l) - α^2 * sin(l) + 2*h*k*cos(l)
                    2*(h*sin(l) - k*cos(l))]

    v⃗ = 1/s^2 * √(μ/p) *[-(sin(l) +α^2*sin(l) - 2*h*k*cos(l) + g -2*f*h*k + α^2*g)
                            -(-cos(l) + α^2*cos(l) + 2*h*k*sin(l) - f +2*g*h*k + α^2*f)
                            2*(h*cos(l) + k*sin(l) + f*h + g*k)]

    vcat(r⃗, v⃗)

    # smp = √(μ/p)

    # # Initialize intermediate values
    # f̂ = zeros(typeof(k), 3)
    # ĝ = zeros(typeof(k), 3)

    # f̂[1] = 1 - k^2 + h^2
    # f̂[2] = 2*k*h
    # f̂[3] = -2k
    # ĝ[1] = 2*k*h
    # ĝ[2] = 1 + k^2 - h^2
    # ĝ[3] = 2h

    # f̂ /= s^2
    # ĝ /= s^2

    # x = r*cos(l)
    # y = r*sin(l)
    # ẋ = -smp * (g + sin(l))
    # ẏ = smp * (f + cos(l))

    # r⃗ = x*f̂ + y*ĝ
    # v⃗ = ẋ*f̂ + ẏ*ĝ
    # vcat(r⃗, v⃗)
end

function Cartesian2MEE(x⃗; μ)
    """
    https://degenerateconic.com/modified-equinoctial-elements.html
    """
    r⃗ = x⃗[1:3]
    v⃗ = x⃗[4:6]
    h⃗ = cross(r⃗, v⃗)

    r = norm(r⃗)
    v = norm(v⃗)
    h_mag = norm(h⃗)

    r̂ = r⃗/r
    ĥ = h⃗/h_mag
    v̂ = (r*v⃗ - r⋅v*r̂)/h_mag
    v̂ = v⃗/v

    p = h_mag*h_mag/μ
    k = ĥ[1]/(1 + ĥ[3])
    h = -ĥ[2]/(1 + ĥ[3])
    s = get_s(;h=h, k=k)

    e⃗ = cross(v⃗, h⃗)/μ - r̂

    # Initialize intermediate values
    f̂ = zeros(typeof(k), 3)
    ĝ = zeros(typeof(k), 3)


    f̂[1] = 1 - k^2 + h^2
    f̂[2] = 2*k*h
    f̂[3] = -2k
    
    ĝ[1] = 2*k*h
    ĝ[2] = 1 + k^2 - h^2
    ĝ[3] = 2h

    f̂ /= s^2
    ĝ /= s^2
 
    f = e⃗ ⋅ f̂
    g = e⃗ ⋅ ĝ
    L = atan( r̂[2] - v̂[1], r̂[1] + v̂[2] )

 
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

    A = [0;0;0;0;0;sqrt(μ*p)*(q/p)^2]
end

function B_equinoctial(MEE; μ)
    """
    See: https://ai.jpl.nasa.gov/public/documents/papers/AAS-22-015-Paper.pdf
    Eq 3
    """
    
    p, f, g, h, k, l = MEE
    q = get_q(f=f, g=g, L=l)
    s = get_s(h=h, k=k)
    
    B = [               0                2*p/q*sqrt(p/μ)                                    0
         sqrt(p/μ)*sin(l)  sqrt(p/μ)*q*((q+1)*cos(l) + f) -sqrt(p/μ)*g/q*(h*sin(l) - k*cos(l))
        -sqrt(p/μ)*cos(l)  sqrt(p/μ)*q*((q+1)*sin(l) + g)  sqrt(p/μ)*f/q*(h*sin(l) - k*cos(l))
                        0                               0             sqrt(p/μ)*s*cos(l)/(2*q)
                        0                               0             sqrt(p/μ)*s*sin(l)/(2*q)
                        0                               0  sqrt(p/μ)*1/q*(h*sin(l) - k*cos(l))]
end


function get_control(MEE; params) # Get control thrust direction and magnitude [0, 1]
    μ = params.μ # problem parameters

    x⃗ = MEE2Cartesian(MEE; μ)
    v⃗ = view(x⃗, 4:6)
    R_inrt2lvlh = DCM_inertial_to_lvlh(x⃗)

    # Just roll with tangential firing to sanity check
    û_LVLH = R_inrt2lvlh*v⃗/norm(v⃗) # thrust along velocity vector (in LVLH frame)
    δ = 1         # Full throttle
    
    û_LVLH, δ 
end

function EOM_MEE!(ẋ, x, p, t)
    """
    See: https://ai.jpl.nasa.gov/public/documents/papers/AAS-22-015-Paper.pdf
    Eq 1
    """
    μ, c, T = p.μ, p.c, p.T_max # problem parameters
    MEE = x[1:6]
    m   = x[7]

    A = A_equinoctial(MEE; μ=μ)
    B = B_equinoctial(MEE; μ=μ)
    û, δ  = get_control(MEE; params=p) # Get control thrust direction and magnitude [0, 1]

    MEE_dot = A + T*δ/m * B * û
    ṁ = -δ*T/(c)

    ẋ[:] = vcat( MEE_dot, ṁ )

end