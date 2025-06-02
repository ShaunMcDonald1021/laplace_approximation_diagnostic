This repo contains all of the code necessary to replicate the results and figures from the manuscript "A probabilistic diagnostic for Laplace approximations: Introduction and experimentation".

Users must first follow the instructions to set up the code from [this repo](https://github.com/tskarvone/fskq) in order for the MATLAB functionality to work.

Other than that, the code is fully self-contained: all data and functions used in the manuscript can be recovered simply by running it, and we have set the seeds such that the results should be identical. The code files here are as follows.

1. `banana.m`: this code generates the results from Section 6.2 of the manuscript.
2. `curse_dim_demo.m`: this generates Figure 8.
3. `diag_calib.m`: a MATLAB function to heuristically calibrate the diagnostic. Users can supply a preliminary grid and candidate values for the free hyperparameters, and see how they affect the results on the chosen calibration function.
4. `fisheries_demo_figs.m`: generate the figures from Section 7.2. Users must run `laplace_diag_demo.R` first to generate the necessary data.
5. `fisheries_demo_timing.m`: compare the timing of the various methods applied to the fisheries data in Section 7. Users must run `GH_stuff_timing.m` and `laplace_diag_demo.R` first.
6. `GH_stuff_timing.m`: timing of the computation of the Gauss-Hermite grids used on the fisheries data. This computation is done separately, so we time it separately.
7. `importace_sampling.R`: anciliary code used to check the validitity of the importance sampling done in Section 7.
8. `L2_minimization_2d.m`: this code allows users to find approximate L2-optimal hyperparameters in 2 dimensions, as described in Section 5.1.
9. `lap_diag.m`: the main function that actually contains the code to output the diagnostic values for a given integrand.
10. `laplace_diag_demo.R`: R code to run the experiments of Section 7 on the fisheries data (which is included in the relevant R package and automatically loaded in the code). Users must first run `GH_stuff_timing.m`.
11. `laplace_diag_functions.R`: helper functions for the above.
12. `laplace_ingredients.m`: uses symbolic math to help calculate the Laplace approximation of a function.
13. `skew_bimodal.m`: this code generates the results from Section 6.1.
14. `visual_calibration_app.mlapp`: this app creates adjustable figures for use in low-dimensional diagnostic calibration, as explored in Section 5.1.
