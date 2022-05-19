%% This script simulates several data sets according to the FSI model and executes the corresponding estimation for each
% Two estimators are considered - the FSI estimator and the Local Frechet
% estimator generalized to p covariates

clear all; close all;
usePar = 0; % should simulations be run in parallel? (1 = yes, 0 = no)
numWk = 6; % if usePar = 1, how many workers should be used.
optRet = 0; % if 0, optimization notes for FSI estimation will not be returned; set to 1 if such notes are desired.

% Adding path to access files in 'manopt' folder
addpath(genpath('manopt'))
% Adding path to access files in 'FMINSEARCHBND' folder for Nelder-Mead optimization
addpath(genpath('FMINSEARCHBND'))


%% Simulation Settings
% These first 4 variables will appear in the filename of the results
nsim = 5;    % number of simulations
n = 50;        % number of observations 
p = 2;         % dimension of index parameter 
noise = 'High'; % High or low noise setting
H = 5; % length of bandwidth that will be created

% Parameters used to initialize the optimization for theta
initType = 'fixed'; % 'random' or 'fixed' - how should the initialization grid be chosen?
nsp = 4; % number of starting points to be used in each angular component, if initType == 'fixed'
nrnd = 4; %number of random starting points to generate, if initType == 'random'

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
b = tmp/norm(tmp); % normalize to have norm 1 (this is theta_0 in the paper)

%% Generate Data

rng(s_dt);      % setting seed for generating data
  
% Generate Uniform(0,1) random variables in matrix n x nsim x p
  % dimensions, they are used as covariates, p gives us the number of
  % covariates used in the model. Hence the data space is (0,1)^p
  
  x = rand(n, nsim, p);  % generate the explanatory variables
  
% generate the index (linear combinations) with given parameters, this
  % gives us the matrix of dimension n x nsim , we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), as required by the function mreg below.  
  
  z = cell2mat(arrayfun(@(j) squeeze(x(:, j, :))*b, 1:nsim, 'UniformOutput', false))/sqrt(p);
  
% Generating the data (x,Y) from fixed index parameter value above for a given p, we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

  % Map the point in (-1,1) on to the circumference of unit sphere, this
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
        for j = (i + 1):n
            normMat(i, j) = norm(squeeze(x(i, k, :) - x(j, k, :)));
            normMat(j, i) = normMat(i, j);
        end
    end

    normSt = sort(normMat, 2);
    bw_min = max(bw_min, max(normSt(:, 4))); bw_max = max(bw_max, max(normSt(:, n)));

end    

% make bandwidth grid that is linear on the log scale
h = exp(linspace(log(bw_min*1.25), log(bw_max), H));

%% Perform Model Fitting

fsiFitAll = cell(1, nsim);
LFpcovFitAll = fsiFitAll;
  
if(usePar == 1)
    
    pool = parpool(numWk);

    parfor i = 1:nsim

        fsiFitAll{i} = get_sphere_fit_FSI(Y{i}, squeeze(x(:, i, :)), h, [], theta_init, optRet);
        LFpcovFitAll{i} = get_sphere_fit_LFpcov(Y{i}, squeeze(x(:, i, :)), h, []);

    end

    delete(pool);

else

    for i = 1:nsim
        
        disp(['Running estimation for dataset ', num2str(i), ' of ', num2str(nsim), ' total simulations'])

        fsiFitAll{i} = get_sphere_fit_FSI(Y{i}, squeeze(x(:, i, :)), h, [], theta_init, optRet);
        LFpcovFitAll{i} = get_sphere_fit_LFpcov(Y{i}, squeeze(x(:, i, :)), h, []);

    end

end

%{
% the MSPE and eta hats corresponding to each simulation are stored as
% matrices instead of struct array for further analyses.
mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
% we find the best bandwidth over all the simulated data. 
mspe_nm_sum = sum(mspe_nm, 2);
h_opt=h(find(mspe_nm_sum == min(mspe_nm_sum(:))));

% mean and standard deviation of the MSPE over simulations at the best bandwidth~
mean(mspe_nm(find(mspe_nm_sum == min(mspe_nm_sum(:))), :));
std(mspe_nm(find(mspe_nm_sum == min(mspe_nm_sum(:))), :));

% to obtain distribution of beta estimates we convert the etas to from
% polar to cartesian coordinates 

eta_opt_s  = zeros(nsim, p-1);
beta_opt_s = zeros(nsim, p);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
    beta_opt_s(i,:)= polar2cart(eta_opt_s(i,:),1);
end


% MSE of the estimator over simulations, their mean and standard deviation
mean(acos(beta_opt_s*b'))
std(acos(beta_opt_s*b'))
 
% computing MSEE from the estimated parameters~
zt = zeros(n, nsim);
Yhat = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
msee = zeros(200, 1);

for i=1:nsim
    disp(i)
    if(p==2)
    zt(:,i) = [x(:,i,1), x(:,i,2)]*beta_opt_s(i,:)';
    
    elseif(p==3)
    zt(:,i) = [x(:,i,1), x(:,i,2), x(:,i,3)]*beta_opt_s(i,:)';
   
    elseif(p==4)
    zt(:,i) = [x(:,i,1), x(:,i,2), x(:,i,3), x(:,i,4)]*beta_opt_s(i,:)';
    
    end
    Yhat{i}= get_sphere_fit_FSI(Y{i}, zt(:,i), h_opt);
    msee(i)= mean((arrayfun(@(g) acos(mreg{i}(:, g)'*Yhat{i}(:,g)),1:n)).^2);
end

% find the mean and standard deviation of MSEE
mean(msee)
std(msee)
%}

fnm = strcat('NM_Sphere_results_n', num2str(n), '_nsim', num2str(nsim), '_p', num2str(p), '_noise', noise, '_', init_Type, '.mat');

save(fnm, 'fsiFitAll', 'LFpcovFitAll', 'n', 'p', 'nsim', 'tau', 's_pm', 's_dt', 'nsp', 'b', 'h_FSI', 'h_LFpcov', 'theta_init', 'usePar', 'numWk')




