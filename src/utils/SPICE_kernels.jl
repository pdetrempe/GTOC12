using Downloads: download

function download_all_kernels()
    const LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
    const SPK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
    const PCK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc"

    ## Only call downloads if you don't already have the SPICE kernels
    # # Download kernels
    download(LSK, "naif0012.tls")
    download(PCK, "gm_de440.tpc")
    download(SPK, "de440.bsp")
end

function furnish_all_kernels()
    # Load leap seconds kernel
    furnsh("naif0012.tls")

    # Load a planetary ephemeris kernel
    furnsh("de440.bsp")

    # Load a planetary properties kernel
    furnsh("gm_de440.tpc")
end 