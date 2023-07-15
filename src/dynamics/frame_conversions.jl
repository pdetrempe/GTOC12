export DCM_inertial_to_rtn

function DCM_inertial_to_rtn(x⃗_inrt)
    r⃗ = view(x⃗_inrt, 1:3)
    v⃗ = view(x⃗_inrt, 4:6)
    h⃗ = cross(r⃗, v⃗)

    r̂ = r⃗/norm(r⃗)
    ĥ = h⃗/norm(h⃗)

    ŷ = cross( ĥ, r̂ )

    R_inrt2rtn = [r̂'; ŷ'; ĥ']
end

function DCM_inertial_to_lvlh(x⃗_inrt)
    r⃗ = view(x⃗_inrt, 1:3)
    v⃗ = view(x⃗_inrt, 4:6)
    h⃗ = cross(r⃗, v⃗)

    r̂ = r⃗/norm(r⃗)
    ĥ = h⃗/norm(h⃗)
    v̂ = v⃗/norm(v⃗)

    x̂ = cross( ĥ, v̂ )

    R_inrt2lvlh = [x̂'; v̂'; ĥ']
end