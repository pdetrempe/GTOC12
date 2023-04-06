# GTOC_1
Repo for goofing around with solving the GTOC_1 problem in Julia

# References
- [Competition](https://sophia.estec.esa.int/gtoc_portal/?page_id=13)
- [Problem Description](https://sophia.estec.esa.int/gtoc_portal/wp-content/uploads/2012/11/ACT-MEM-MAD-GTOC1-The-Problem_V4.pdf)

# Gameplan
0. Plotting
   1. Show some trajectory in the Sun-centered J2000 frame
1. Execute feasible solution
   1. Dynamics models
      1. Ephemeris
      2. Spacecraft
         1. Sun gravity
         2. Flybys
            1. Minimum Perilune constraint
         3. Actuator model
      3. Asteroid
         1. Closed-form Keplerian orbit
   2. Evaluate Objective Function
   3. Use the two above to recreate the winning solution
2. Converge perturbed winning solution with shooting methods
3. Add local optimizer to converged method
4. Wrap shooting+local optimizer in global hyperparameter optimizer
