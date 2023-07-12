
export download_all_kernels, furnish_all_kernels, get_planet_state

function download_all_kernels()
    LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
    SPK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
    PCK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc"

    ## Only call downloads if you don't already have the SPICE kernels
    # # Download kernels
    if !isdir("deps")
        mkdir("deps")
    end
    
    download(LSK, "./deps/naif0012.tls")
    download(PCK, "./deps/gm_de440.tpc")
    download(SPK, "./deps/de440.bsp")
end

function furnish_all_kernels()
    # Check if kernels exist first and download them if need be
    if !isdir("deps") || !isfile("./deps/naif0012.tls")
        print("Downloading necessary SPICE kernels")
        download_all_kernels()
    end

    # Load leap seconds kernel
    furnsh("./deps/naif0012.tls")

    # Load a planetary ephemeris kernel
    furnsh("./deps/de440.bsp")

    # Load a planetary properties kernel
    furnsh("./deps/gm_de440.tpc")
end 

# Wrap SPICE calls
function get_planet_state(
    planet::Union{String,Int}, 
    ET::AbstractFloat; 
    ref_frame::String=default_ref_frame, 
    observer::Union{String,Int}=default_CB_idx
    )
    _planet = if planet isa AbstractString
        bodn2c(planet)
    elseif planet isa Integer
        planet
    end
    _observer = if observer isa AbstractString
        bodn2c(observer)
    elseif observer isa Integer
        observer
    end
    r_planet = spkgeo(_planet, ET, ref_frame, _observer)[1]
    r_planet *= 1000 # km → m
end