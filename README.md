# GTOC12
Repo for (trying to) solve the GTOC12 problem in Julia

- [Competition](https://gtoc12.tsinghua.edu.cn/competition)
- [Problem Description](https://gtoc12.tsinghua.edu.cn/competition/theProblem)

# References
- [NASA JPL NAIF](https://naif.jpl.nasa.gov/naif/)
   - [NAIF Leapseconds Kernel (LSK)](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/)
   - [NAIF Planetary Constants Kernel (PCK)](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/)
      - [NAIF PCK Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html)
   - [NAIF Spacecraft Planet Kernel (SPK)](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/)
      - [NAIF SPK Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [JuliaMono font](https://juliamono.netlify.app/)
## Trajectory Optimization
- [ESA PyKEP](https://github.com/esa/pykep)
- [shooting methods reference](https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2013_Pavlak.pdf)
- [Low-Thrust Trajectory Optimization for the Solar System Pony Express](https://ai.jpl.nasa.gov/public/documents/papers/AAS-22-015-Paper.pdf)
- [Multi-Gravity Assist Low-Thrust Trajectory Optimization](https://www.researchgate.net/profile/Martin-Ozimek/publication/330533150_The_Low-thrust_Interplanetary_EXplorer_A_Medium-Fidelity_Algorithm_for_Multi-Gravity_Assist_Low-Thrust_Trajectory_Optimization/links/5c46900e92851c22a3870254/The-Low-thrust-Interplanetary-EXplorer-A-Medium-Fidelity-Algorithm-for-Multi-Gravity-Assist-Low-Thrust-Trajectory-Optimization.pdf)
## Modified Equinoctial Orbit Elements
- [Modified Equinoctial Elements](https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf)
- [Fortran Astrodynamics Toolkit, Modified Equinoctial Module](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit/blob/master/src/modified_equinoctial_module.f90)
- [Degenerate Conic | Modified Equinoctial Elements](https://degenerateconic.com/modified-equinoctial-elements.html)
- [*A Set of Modified Equinoctial Elements*, Walker, 1985](http://cdsads.u-strasbg.fr/pdf/1985CeMec..36..409W)
   - [*Errata - A Set of Modified Equinoctial Elements*, Walker, 1986](http://cdsads.u-strasbg.fr/pdf/1986CeMec..38..391W)
- [*On the Equinoctial Orbit Elements*, Broucke, 1972](https://adsabs.harvard.edu/full/1972CeMec...5..303B)

# Gameplan
0. Plotting
   - [x] Show some trajectory in the Sun-centered J2000 frame
1. Execute feasible solution
   1. Dynamics models
      - [x] Ephemeris
      - [ ] Good ol Kepler
      2. Spacecraft
         - [x] Sun gravity
         - [ ] Flybys
            - [ ] Minimum Perilune constraint
         - [x] Launch vehicle
         - [x] Actuator model
      - [x] Asteroid
         - [x] Closed-form Keplerian orbit
         - [ ] Test Keplerian time/mean anomaly conversion
   2. Evaluate Objective Function
   3. Use the two above to recreate the winning solution
2. Converge trajectory to an asteroid with shooting methods
3. Add local optimizer to converged method
4. Wrap shooting+local optimizer in global hyperparameter optimizer
5. Start doing tree searches for global optima
