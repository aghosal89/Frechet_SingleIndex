The manifold optimization toolbox 'Manopt' is available to download from the link:
https://www.manopt.org/download.html

After the folder 'manopt' is extracted from the zip file 'Manopt_6.0', a path is generated in 'updated_codes_frechet.m' and 'run_sphere_sim_neldermead.m' to work with various functions necessary for computation. 

Running the codes in the 'updated_codes_frechet.m' requires creating a directory folder in Matlab which contains 'manopt' folder above and the following functions:

1) add_noise.m
2) cart2polar.m
3) get_cost.m
4) get_egrad.m
5) get_ehess.m
6) get_sphere_fit_SIM2.m
7) polar2cart.m

Running the codes in 'run_sphere_sim_neldermead.m' requires adding the function nm_cost.m 

to the folder above 
