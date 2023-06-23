using LinearAlgebra

export DCM_inertial_to_lvlh

function DCM_inertial_to_lvlh(x⃗_inrt)
    r⃗ = view(x⃗_inrt, 1:3)
    v⃗ = view(x⃗_inrt, 4:6)
    h⃗ = cross(r⃗, v⃗)

    r̂ = r⃗/norm(r⃗)
    ĥ = h⃗/norm(h⃗)

    ŷ = cross( ĥ, r̂ )

    R_inrt2lvlh = [r̂'; ŷ'; ĥ']
end