The manifold optimization toolbox 'Manopt_7.0' is available to download from the link:
https://www.manopt.org/download.html

After the folder 'manopt 2' is extracted from the zip file 'Manopt_7.0', a path is generated in the script 'run_sphere_sim_neldermead.m' to work with various functions necessary for computation. 

For implementing Neldermead optimization download the zip file 'FMINSEARCHBND.zip' extract the folder 'FMINSEARCHBND' from the following link:
https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon , to be used in the simulation, a path is generated in 'run_sphere_sim_neldermead.m' in same way as above.

Running the codes in the 'run_sphere_sim_neldermead.m' requires creating a directory folder in Matlab which contains 'manopt 2' folder above and the following functions:

- add_noise.m
- cart2polar.m
- get_cost.m
- get_egrad.m
- get_ehess.m
- get_sphere_fit_FSI.m
- polar2cart.m
- nm_cost.m
- fminsearchbnd.m

Once the optimization steps are completed for various simulation settings, the estimates of the index parameter and the choice of bandwidth are computed using codes in:

- sphere_computations.m

And the necessary plots are generated using:

- sphere_plots.m

