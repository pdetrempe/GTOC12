using SPICE

export download_all_kernels, furnish_all_kernels

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
function get_planet_state(planet, ET)
    # Use SPICE to get position
    ref = "ECLIPJ2000"
    abcorr = "NONE"
    obs = GTOC12.default_CB_str

    # TODO: Use SPICE Int IDs instead of string
    r_planet = spkezr(uppercase(planet), ET, ref, abcorr, uppercase(obs))[1]
    r_planet *= 1000 # km → m
end