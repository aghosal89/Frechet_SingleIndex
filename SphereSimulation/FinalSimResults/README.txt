This folder contains 18 matlab data (.mat) files corresponding to 18 different
simulation settings: 3 sample sizes (n = 50, 100, 200), three covariate dimension
values (p = 2, 5, 10) and two noise levels (Low and High).  For example, the file

Sphere_results_n50_nsim200_p2_noiseHigh_Final.mat

stores the results for sample size n = 50 and dimension p = 2 under low noise. All
simulation settings include 200 independent simulation runs.

Each file contains the following variables that can be used to visualize results
and compute summary statistics.

n - sample size (integer)
p - covariate dimension (integer)
theta0 - true coefficient parameter used to generate the data under the FSI model;
         (vector of length p)
h - vector of bandwidths used in estimation (vector of length H, H depending on the setting)
tau - noise standard deviation (numeric)
x - covariate values for all data sets (numeric array of size n-by-200-by-p); hence,
    x(:, m, :) represents the n-by-p covariate matrix for the m-th data set
mreg - true conditional Frechet means, i.e. the values of the true regression curve
       evaluated at the covariate values in x, on the sphere (cell array of length 200);
       mreg{m} is a 3-by-n matrix with i-th column giving the value of the regression
       function evaluated at x(i, m, :)
Y - response values corresponding to x (cell array of length 200); Y{m} is a 3-by-n
    matrix with i-th column giving the response value for the covariate vector x(i, m, :)
LFpcovFitAll - results obtained by multivariate local Frechet regression (length 200 cell array);
               LFpcovFitAll{m} will be a 3-by-n-by-H array of fitted values, with
               the third dimension indexing the different bandwidths.  For further
               detail, see the file get_sphere_fit_LFpcov.m
fsiFitAll - results obtained by fitting the FSI model (length 200 cell array), including
            estimates of theta0 as well as the fitted values.  For more detail, see
            the documentation in get_sphere_fit_FSI.m
