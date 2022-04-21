1) The manifold optimization toolbox for matlab used for our analyses, 'Manopt_6.0' is available to download from the following link:

https://github.com/NicolasBoumal/manopt-web/blob/master/downloads/Manopt_6.0.zip

Unzip 'Manopt_6.0', which will give folder 'manopt'. Move folder 'manopt' into the directory from which you are running Matlab.

2) For implementing Neldermead optimization download the zip file 'FMINSEARCHBND.zip' from the following link:
https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon

Extract the folder 'FMINSEARCHBND' by unzipping 'FMINSEARCHBND.zip' then move it to the same directory as the 'manopt'. 

3) Open Matlab and set the working directory to the folder containing the folders 'manopt' and 'FMINSEARCHBND'. 
a path is generated in the scripts:

4) download from https://github.com/aghosal89/Frechet_SingleIndex/tree/master/Sphere%20Simulation the scripts:

  - run_sphere_SIM_neldermead.m
  - run_sphere_SIM_grid.m 
  - add_noise.m
  - cart2polar.m
  - get_cost.m
  - get_egrad.m
  - get_ehess.m
  - get_sphere_fit_FSI.m
  - polar2cart.m
  - nm_cost.m

5)a) run 'run_sphere_SIM_neldermead.m' 
To simulate the data and produce numbers to fill figures... and tables... in the paper. 

5)b) and 'run_sphere_SIM_grid.m'
Running the codes in the 'run_sphere_SIM_neldermead.m' requires the following files in the directory folder which also contains folders 'manopt' and 
'FMINSEARCHBND' mentioned above:

To run the codes in the script 'run_sphere_SIM_grid.m' the above functions are enough as kept in the same working directory.  
