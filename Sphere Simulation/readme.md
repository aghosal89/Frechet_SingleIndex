The manifold optimization toolbox for matlab used for our analyses, 'Manopt_6.0' is available to download from the following link:

https://github.com/NicolasBoumal/manopt-web/blob/master/downloads/Manopt_6.0.zip

The folder 'manopt' is extracted from the zip file 'Manopt_6.0' and stored in the working directory.

For implementing Neldermead optimization download the zip file 'FMINSEARCHBND.zip' extract the folder 'FMINSEARCHBND' from the following link:
https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon , 

a path is generated in the scripts:

1)run_sphere_SIM_neldermead.m
2)run_sphere_SIM_grid.m 

Running the codes in the 'run_sphere_SIM_neldermead.m' requires the following files in the directory folder which also contains folders 'manopt' and 
'FMINSEARCHBND' mentioned above:

- add_noise.m
- cart2polar.m
- get_cost.m
- get_egrad.m
- get_ehess.m
- get_sphere_fit_FSI.m
- polar2cart.m
- nm_cost.m

To run the codes in the script 'run_sphere_SIM_grid.m' the above functions are enough as kept in the same working directory.  
