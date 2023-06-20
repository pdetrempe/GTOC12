export download_all_kernels, furnish_all_kernels

function download_all_kernels()
    LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
    SPK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
    PCK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc"

    ## Only call downloads if you don't already have the SPICE kernels
    # # Download kernels
    download(LSK, "./deps/naif0012.tls")
    download(PCK, "./deps/gm_de440.tpc")
    download(SPK, "./deps/de440.bsp")
end

function furnish_all_kernels()
    # Load leap seconds kernel
    furnsh("./deps/naif0012.tls")

    # Load a planetary ephemeris kernel
    furnsh("./deps/de440.bsp")

    # Load a planetary properties kernel
    furnsh("./deps/gm_de440.tpc")
end 