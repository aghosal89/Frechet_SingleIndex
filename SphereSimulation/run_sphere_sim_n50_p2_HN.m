%% This script simulates several data sets according to the FSI model and executes the corresponding estimation for each
% Two estimators are considered - the FSI estimator and the Local Frechet
% estimator generalized to p covariates

clear all; close all;
usePar = 1; % should simulations be run in parallel? (1 = yes, 0 = no)
numWk = 8; % if usePar = 1, how many workers should be used.

% Adding path to access files in 'manopt' folder
addpath(genpath('manopt'))
% Adding path to access files in 'FMINSEARCHBND' folder for Nelder-Mead optimization
addpath(genpath('FMINSEARCHBND'))


%% Simulation Settings
% These first 4 variables will appear in the filename of the results
nsim = 50;    % number of simulations
n = 50;        % number of observations 
p = 2;         % dimension of index parameter 
noise = 'High'; % High or low noise setting
H = 10; % length of bandwidth that will be created

% Parameters used to initialize the optimization for theta
initType = 'random'; % 'random' or 'fixed' - how should the initialization grid be chosen?
nsp = 3; % number of starting points to be used in each angular component, if initType == 'fixed'
nrnd = 4; % number of random starting points to generate, if initType == 'random'


% Optimization options -  If optType = 'TR', runs trust region (fmincon)
% with analytic gradient/Hessian of approximate 
% cost function using local linear instead of local Frechet estimates to
% increase speed.  Otherwise, Nelder-Mead (fminsearchbnd) is used
optType = 'TR'; 
if(strcmp(optType, 'TR'))
    optns = struct('Algorithm', 'TR', 'Ret', 1, 'useAltGrad', 1, 'useAltHess', 1, 'maxIter', 100); % TR with gradient and Hessian approximation and max of 100 iterations per starting value
else
    optns = struct('Algorithm', 'NM', 'Ret', 1, 'maxIter', 100); % NM with max of 100 iterations per starting value
end

% Define random seeds
s_pm = 19878;  % seed for generating true theta parameter (should be same across all sims/batches sharing same p)
s_dt = 285720;  % seed for generating random data (different values for different sims/batches)

%% Population parameters for data generation

% Variance for data simulation
if strcmp(noise, 'High')
    tau = sqrt(0.6);    % Noise
else
    tau = sqrt(0.2);
end

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true theta for the number of covariates p:
rng(s_pm);      % setting seed for generating index parameter

tmp = 2*rand(p, 1) - 1; % p uniform rvs between -1 and 1
theta0 = tmp/norm(tmp); % normalize to have norm 1

%% Generate Data

rng(s_dt);      % setting seed for generating data
  
% Generate Uniform(0,1) random variables in matrix n x nsim x p
  % dimensions, they are used as covariates, p gives us the number of
  % covariates used in the model. Hence the data space is (0,1)^p
  
  x = rand(n, nsim, p);  % generate the explanatory variables
  
% generate the index (linear combinations) with given parameters, this
  % gives us the matrix of dimension n x nsim , we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), as required by the function mreg below.  
  
  z = cell2mat(arrayfun(@(j) squeeze(x(:, j, :))*theta0, 1:nsim, 'UniformOutput', false))/sqrt(p);
  
% Generating the data (x,Y) from fixed index parameter value above for a given p, we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

  % Map the point in (-1,1) onto the unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 

  % Generate N(0, (tau)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  % These represent random errors in the tangent space 
  v1 = tau*randn(n, nsim); 
  v2 = tau*randn(n, nsim); 

 
  % create the object for storing the simulated response data based on
  % our assumptions on theta and covariates: 
  
  Y = arrayfun(@(i) cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j), [v1(j, i) v2(j, i)]), 1:n, 'UniformOutput', false)), 1:nsim, 'UniformOutput', false);
  
%% Create grid of starting values for theta and grids of bandwidths

% Creating the grid of polar coordinates as starting points for
% p-dimensional hypersphere, stored in theta_init

if(strcmp(initType, 'fixed'))
    theta_init = thetaGrid(nsp, p - 1);
else    
    theta_init = nrnd;
end

% choosing a range of bandwidths so that, for all simulated data sets and 
% for all i = 1,...n, the i-th predictor (a p-vector) will have at least 3 
% neighbors within a distance equal to the smallest bandwidth.  Note that
% this also implies that, for any p-vector theta on the unit sphere, the
% projected covariates (in the direction theta) will also have at least 3
% neighbors within this distance.  The largest bandwidth is set to be the
% largest difference, across simulations, between any two predictors within
% a simulation.

bw_min = 0;
bw_max = 0;

for k = 1:nsim

    normMat = zeros(n, n);
    for i = 1:(n - 1)
        for j = (n + 1):n
            normMat(i, j) = norm(squeeze(x(i, :, k) - x(j, :, k)));
            normMat(j, i) = normMat(i, j);
        end
    end

    normSt = sort(normMat, 2);
    bw_min = max(bw_min, max(normSt(:, 4))); hMax = max(bw_max, max(normSt(:, n)));

end    

% make bandwidth grid that is linear on the log scale
h = exp(linspace(log(bw_min*1.05), log(bw_max), H));

%% Perform Model Fitting

fsiFitAll = cell(1, nsim);
LFpcovFitAll = fsiFitAll;
  
if(usePar == 1)
    
    pool = parpool(numWk);

    parfor i = 1:nsim

        fsiFitAll{i} = get_sphere_fit_FSI(Y{i}, squeeze(x(:, i, :)), h_FSI, [], theta_init, optns);
        LFpcovFitAll{i} = get_sphere_fit_LFpcov(Y{i}, squeeze(x(:, i, :)), h_LFpcov, []);

    end

    delete(pool);

else

    for i = 1:nsim

        fsiFitAll{i} = get_sphere_fit_FSI(Y{i}, squeeze(x(:, i, :)), h_FSI, [], theta_init, optns);
        LFpcovFitAll{i} = get_sphere_fit_pcov(Y{i}, squeeze(x(:, i, :)), h_LFpcov, []);

    end

end

fnm = strcat('NM_Sphere_results_n', num2str(n), '_nsim', num2str(nsim), '_p', num2str(p), '_noise', noise, '_', init_Type, '.mat');

save(fnm, 'fsiFitAll', 'LFpcovFitAll', 'n', 'p', 'nsim', 'tau', 's_pm', 's_dt', 'nsp', 'theta0', 'h_FSI', 'h_LFpcov', 'theta_init', 'usePar', 'numWk')






