%% This script simulates several data sets according to the FSI model and executes the corresponding estimation for each
% Two estimators are considered - the FSI estimator and the Local Frechet
% estimator generalized to p covariates

usePar = 1; % should simulations be run in parallel? (1 = yes, 0 = no)
numWk = 20; % if usePar = 1, how many workers should be used.

% Adding path to access files in 'manopt' folder
addpath(genpath('manopt'))
% Adding path to access files in 'FMINSEARCHBND' folder for Nelder-Mead optimization
addpath(genpath('FMINSEARCHBND'))


%% Simulation Settings
% These first 4 variables will appear in the filename of the results
nsim = 20;    % number of simulations
n = 200;        % number of observations
p = 10;         % dimension of index parameter
noise = 'High'; % High or low noise setting
h = 0.15:0.1:0.45; % grid of bandwidths

% Define random seeds
s_pm = 54639;  % seed for generating true theta parameter (should be same across all sims/batches sharing same p)
s_dt = 239446;  % seed for generating random data (different values for different sims/batches)

%% Population parameters for data generation

% Variance for data simulation
if strcmp(noise, 'High')
    tau = sqrt(0.8);    % Noise
else
    tau = sqrt(0.4);
end

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true theta for the number of covariates p:
rng(s_pm);      % setting seed for generating index parameter

tmp = 2*rand(p, 1) - 1; % p uniform rvs between -1 and 1
if(tmp(1) < 0) % make first component positive
  tmp = -tmp;
end
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

%% Perform Model Fitting

% Optimization options - See get_sphere_fit_FSI.m for details
% If optType is 'IP' or 'SQP', runs fmincon in final optimization step.
% Otherwise, Nelder-Mead (fminsearchbnd) is used
optType = 'SQP';
q = 100; % # of random starting points
eta_init = pi*rand(p - 1, q) - pi/2;
thetaInit = cell2mat(arrayfun(@(j) polar2cart(eta_init(:, j)), 1:q, 'UniformOutput', false));

optns = struct('Ret', 1, 'maxIter', 50, 'thetaInit', thetaInit, 'thetaChoose', 5, 'Algorithm', optType); % max of 70 iterations per starting value

fsiFitAll = cell(1, nsim);
LFpcovFitAll = fsiFitAll;
bad = zeros(1, nsim);

tmFSI = zeros(1, nsim);
tmLF = zeros(1, nsim);

if(usePar == 1)

    pool = parpool(numWk);

    parfor i = 1:nsim
        
        tic;
        try
            fsiFitAll{i} = get_sphere_fit_FSI(Y{i}, squeeze(x(:, i, :))/sqrt(p), h, [], optns);
        catch errTmp
            disp(['Error with simulated data set ' num2str(i)])
            bad(i) = 1;
        end
        tmFSI(i) = toc;
        
        tic;
        LFpcovFitAll{i} = get_sphere_fit_LFpcov(Y{i}, squeeze(x(:, i, :))/sqrt(p), h, []);
        tmLF(i) = toc;
        
        disp(['Finished simulation ' num2str(i)])

    end

    delete(pool);

else

    for i = 1:nsim
        
        tic;
        try
            fsiFitAll{i} = get_sphere_fit_FSI(Y{i}, squeeze(x(:, i, :))/sqrt(p), h, [], optns);
        catch errTmp
            disp(['Error with simulated data set ' num2str(i)])
            bad(i) = 1;
        end
        tmFSI(i) = toc;
        
        tic;
        LFpcovFitAll{i} = get_sphere_fit_LFpcov(Y{i}, squeeze(x(:, i, :))/sqrt(p), h, []);
        tmLF(i) = toc;

        disp(['Finished simulation ' num2str(i)])

    end

end

fnm = strcat('Sphere_results_n', num2str(n), '_nsim', num2str(nsim), '_p', num2str(p), '_noise', noise, 'batch1.mat');

save(fnm, 'x', 'Y', 'mreg', 'fsiFitAll', 'LFpcovFitAll', 'n', 'p', 'nsim', 'tau', 's_pm', 's_dt', 'theta0', 'h', 'optns', 'usePar', 'numWk', 'bad', 'tmFSI', 'tmLF')
