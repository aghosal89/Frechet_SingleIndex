%% This script simulates several data sets according to the FSI model and executes the corresponding estimation for each

clear all; close all;
usePar = 1; % should simulations be run in parallel? (1 = yes, 0 = no)
numWk = 6; % if usePar = 1, how many workers should be used.

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

s_pm = 19878;  % seed for generating theta parameter (should be same across all sims/batches sharing same p)
s_dt = 285720;  % seed for generating random data (different values for different sims/batches)
nsp = 3; % number of starting points to be used in each angular component

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

% Generate N(0, (tau)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
% These represent random errors in the tangent space 
  v1 = tau*randn(n, nsim); 
  v2 = tau*randn(n, nsim); 
  
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
 
 
  % create the object for storing the simulated response data based on
  % our assumptions on theta and covariates: 
  
  Y = arrayfun(@(i) cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j), [v1(j, i) v2(j, i)]), 1:n, 'UniformOutput', false)), 1:nsim, 'UniformOutput', false);
  
%% Create grid of starting values for theta and grids of bandwidths

% Creating the grid of polar coordinates as starting points for
% p-dimensional hypersphere, stored in theta_init

theta_init = thetaGrid(nsp, p - 1);

% choosing a range of bandwidths for FSI
bw_min = 0;
bw_max = 0;

for k = 1:size(theta_init, 2)
    for i = 1:nsim
        
        temp = squeeze(x(:, i, :))*theta_init(:, k);
        
        tempd = abs(repmat(temp, 1, length(temp)) - repmat(temp', length(temp), 1));
        tempc = sort(tempd, 2); % nxn 
        bw_min= max(bw_min, max(tempc(:, 4)));
        bw_max= max(bw_max, max(tempc(:, n)));
    end
end    

h_FSI = exp(linspace(log(bw_min*1.25), log(bw_max), 10));

% choosing a range of bandwidths for local Frechet with p covariates
bw_min = 0;
bw_max = 0;

for k = 1:p
    for i = 1:nsim
        
        temp = squeeze(x(:, i, k));
        
        tempd = abs(repmat(temp, 1, length(temp)) - repmat(temp', length(temp), 1));
        tempc = sort(tempd, 2); % nxn 
        bw_min= max(bw_min, max(tempc(:, 4)));
        bw_max= max(bw_max, max(tempc(:, n)));
    end
end    

h_LFpcov = exp(linspace(log(bw_min*1.25), log(bw_max), 10));



%% Perform Model Fitting

fsiFitAll = cell(1, nsim);
LFpcovFitAll = fsiFitAll;
  
if(usePar == 1)
    
    pool = parpool(numWk);

    parfor i = 1:nsim

        fsiFitAll{i} = get_sphere_fit_FSI(Y{i}, squeeze(x(:, i, :)), [], h_FSI, theta_init);
        LFpcovFitAll{i} = get_sphere_fit_pcov(Y{i}, squeeze(x(:, i, :)), [], h_LFpcov);

    end

    delete(pool)

else

    for i = 1:nsim

        fsiFitAll{i} = get_sphere_fit_FSI(Y{i}, squeeze(x(:, i, :)), [], h_FSI, theta_init);
        LFpcovFitAll{i} = get_sphere_fit_pcov(Y{i}, squeeze(x(:, i, :)), [], h_LFpcov);

    end

end
        
%delete(pool) % shutting down parallel pool 
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

fnm = strcat('NM_Sphere_results_n', num2str(n), '_nsim', num2str(nsim), '_p', num2str(p), '_noise', noise, '.mat');

save(fnm, 'fsiFitAll', 'LFpcovFitAll', 'n', 'p', 'nsim', 'tau', 's_pm', 's_dt', 'nsp', 'b', 'h_FSI', 'h_LFpcov', 'theta_init', 'usePar', 'numWk')
%save NM_Sphere_results_n200_nsim200_p2_HN.mat *







