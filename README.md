# GTOC_1
Repo for goofing around with solving the GTOC_1 problem in Julia

# References
- [Competition](https://sophia.estec.esa.int/gtoc_portal/?page_id=13)
- [Problem Description](https://sophia.estec.esa.int/gtoc_portal/wp-content/uploads/2012/11/ACT-MEM-MAD-GTOC1-The-Problem_V4.pdf)

# Gameplan
0. Plotting
   - [x] Show some trajectory in the Sun-centered J2000 frame
1. Execute feasible solution
   1. Dynamics models
      - [x] Ephemeris
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
2. Converge perturbed winning solution with shooting methods
3. Add local optimizer to converged method
4. Wrap shooting+local optimizer in global hyperparameter optimizer
