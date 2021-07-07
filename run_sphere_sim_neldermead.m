%% Computation of MSPE/MSEE for various simulation scenarios using Nelder-Mead optimization

%% Sphere simulations n=50, p=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Adding path to folders
addpath(genpath('manopt'))
addpath(genpath('FMINSEARCHBND'))

% For beta in circumference of a unit circle;

nsim = 10;  % number of simulations
n=100;       % number of observations 
p=2;         % dimension of theta vector 
s=19878;     % random number generators
rng(s);      % setting seed for random number generator

% Variances for normal data simulation
tau1=0.2;   % Low Noise 
tau2=0.6;   % High Noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true theta
b =[cos(-0.704), sin(-0.704)];

% Creating the starting points for eta on hypersphere p=2
start = linspace(-1.57+(pi/15),1.57-(pi/15),5);         

h1 = exp(log(0.35):0.02:log(0.6));     % bandwidth sequence for low noise data
h2 = exp(log(0.45):0.02:log(0.95));    % bandwidth sequence for high noise data

pm = struct('s', s, 'nsim', nsim, 'n', n, 'h1', h1, 'h2', h2, 'tau1', tau1, 'tau2', tau2);

eta_opt_mspe_hn = arrayfun(@(l) zeros(length(h2), p-1, size(start,1)), 1:nsim, 'UniformOutput', false);
eta_opt_mspe_ln = arrayfun(@(l) zeros(length(h1), p-1, size(start,1)), 1:nsim, 'UniformOutput', false);

% Generating the data from fixed \theta_0 value, here we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

  % Generate N(0, (0.2)^2) random numbers in matrix n x nsim in 2 sets u1 and u2 respectively
  u1 = tau1*randn(n, nsim); 
  u2 = tau1*randn(n, nsim); 
  
  % Generate N(0, (0.35)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  v1 = tau2*randn(n, nsim); 
  v2 = tau2*randn(n, nsim); 
  
  % Generate Uniform(0,1) random variables in matrix n x nsim x p
  % dimensions, they are used as covariates, p gives us the number of
  % covariates used in the model. Hence the data space is (0,1)^p
  
  x = rand(n, nsim, p);  
  
  % generate the index (linear combinations) with given parameters, this
  % gives us the matrix of dimension n x nsim , we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), now this z can be fed into the Local frechet regression.
  
  z= (b(1)*x(:,:,1) + b(2)*x(:,:,2))./sqrt(p); 
  
  % Map the point in (-1,1) on to the circumference of unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg1 = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
  % defining a set of functions pre-defined on a unit sphere of p-dimensions
  M = spherefactory(3); % define structures fr, ops to store the pre-defined functions on a sphere
  fr = struct(); 
  ops = struct();
  fr.M = M; 
  ops.verbosity = 0;
 
  % create the variables for storing the simulated response data based on
  % our assumptions on theta and covariates: creates a list of nsim number
  % of matrices of size 3 x n, for each simulation we have 
  
  Y_true1 = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  Y_true2 = Y_true1; 
  
  % for each simulation, generate response data by first transporting the single index
  % onto the unit sphere, then map that point to the tangent space, generates a single point by adding noise created by
  % normal samples in u1, u2.
  
  for i=1:nsim
     Y_true1{i} = cell2mat(arrayfun(@(j) add_noise(mreg1{i}(:, j),[u1(j, i) u2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
  for i=1:nsim
     Y_true2{i} = cell2mat(arrayfun(@(j) add_noise(mreg1{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
% Creating the array with necessary dimensions to store MSPE in both high and low noises respectively 

mspe_ln = arrayfun(@(i) zeros(length(h1), size(start,1)), 1:nsim, 'UniformOutput', false);
mspe_hn = arrayfun(@(i) zeros(length(h2), size(start,1)), 1:nsim, 'UniformOutput', false);

tic % starting the time record for estimation-optimization procedure
pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim % starting the parallel for loop
%for i = 1:nsim
disp(i)
%i=1;
Y_temp1 = Y_true1{i};   % creating temporary response variable to be input for cost function
Y_temp2 = Y_true2{i};   % creating temporary response variable to be input for cost function
xtemp = [x(:,i,1), x(:,i,2)];     % getting the covariates for the i-th simulation
              
for j=1:length(start)
   disp(j)
   for l = 1:length(h1) % for every bandwidth in the low noise
     
       cost_fun = @(eta) nm_cost(Y_temp1, eta, h1(l), xtemp); % define cost function 
       
       % running the optimization for starting point 
       [eta_opt1, fval] = fminsearchbnd(cost_fun, start(j), -1.57, 1.57);
       mspe_ln{i}(l,j) = fval;         % storing minimum cost function
       eta_opt_mspe_ln{i}(l,:,j)= eta_opt1;   % storing local minima
       
    end
    
    for l = 1:length(h2) % bandwidth range for high noise scenario
    
       cost_fun = @(eta) nm_cost(Y_temp2, eta, h2(l), xtemp); % define cost function 
    
       % running the optimization for starting point 
       [eta_opt2, fval] = fminsearchbnd(cost_fun, start(j), -1.57, 1.57);
       mspe_hn{i}(l, j) = fval;          % storing minimum cost function
       eta_opt_mspe_hn{i}(l,:,j) = eta_opt2;   % storing local minima
  
    end
  
end

end

delete(pool) % shutting down parallel pool 
rt=toc; % storing the time elapsed. 

save Sphere_results_n100_nsim10_p2_nm.mat mreg1 * mspe_hn * mspe_ln * eta_opt_mspe_hn * eta_opt_mspe_ln * h1 * h2

%%%%% Code to find the optimum bandwidth~

% start =  -0.83;
% for high noise
mspe_nm_hn = zeros(nsim, length(h2), size(start,1));
% for low noise
mspe_nm_ln = zeros(nsim, length(h1), size(start,1));

for i=1:nsim
    for k=1:size(start,1)
        for j=1:length(h2)
            mspe_nm_hn(i,j,k)= mspe_hn{i}(j,k);
        end
        for j=1:length(h1)
            mspe_nm_ln(i,j,k)= mspe_ln{i}(j,k);
        end
    end
end

    
% to find the average and standard deviation of the MSPE
% high noise~
mean(mspe_nm_hn(:,:,1))

% distribution of beta estimates 
% high noise
beta_opt_mspe_hn_mat = zeros(nsim, p);
for i=1:nsim
    beta_opt_mspe_hn_mat(i,:)= polar2cart(eta_opt_mspe_hn{i},1);
end

% MSE  of theta estimator
mean(acos(beta_opt_mspe_hn_mat*b'))
std(acos(beta_opt_mspe_hn_mat*b'))

% low noise
beta_opt_mspe_ln_mat = zeros(nsim, p);
for i=1:nsim
    beta_opt_mspe_ln_mat(i,:)= polar2cart(eta_opt_mspe_ln{i},1);
end

% MSE  of theta estimator
mean(acos(beta_opt_mspe_ln_mat*b'))
std(acos(beta_opt_mspe_ln_mat*b'))



















%% Sphere simulations n=50, p=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Adding path to folders
addpath(genpath('manopt'))
addpath(genpath('FMINSEARCHBND'))

% For beta in circumference of a unit circle;

nsim = 10;  % number of simulations
n=50;       % number of observations 
p=3;         % dimension of theta vector 
s=19878;     % random number generators
rng(s);      % setting seed for random number generator

%Variances for normal data simulation
tau1=0.2;   % Low noise 
tau2=0.35;  % High Noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true beta parameter
b=[0.386, 0.181, -sqrt(1-0.386^2-0.181^2)];

% Creating the grid for theta on hypersphere p=3, using polar coordinates
%phi = -pi+(pi/5):(pi/1.2): pi-(pi/5); %Choose the grid for eta 
%[S1 S2] = ndgrid(phi', phi');
%start = [S1(:), S2(:)]; % creating the grid

% plot check
%plot(start(:,1), start(:,2), 'r.');
%ylim = [-pi, pi];
%xlim = [-pi, pi];
%size(start)

h1 = exp(log(0.45):0.02:log(0.65));     % bandwidth sequence for low noise data
h2 = exp(log(0.45):0.02:log(0.65));    % bandwidth sequence for high noise data

pm = struct('s', s, 'nsim', nsim, 'n', n, 'h1', h1, 'h2', h2, 'tau1', tau1, 'tau2', tau2);

eta_opt_mspe_hn = arrayfun(@(l) zeros(length(h2),p-1, size(start,1)), 1:nsim, 'UniformOutput', false);
eta_opt_mspe_ln = arrayfun(@(l) zeros(length(h1),p-1, size(start,1)), 1:nsim, 'UniformOutput', false);

% Generating the data from fixed \theta_0 value, here we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

  % Generate N(0, (0.2)^2) random numbers in matrix n x nsim in 2 sets u1 and u2 respectively
  u1 = tau1*randn(n, nsim); 
  u2 = tau1*randn(n, nsim); 
  
  % Generate N(0, (0.35)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  v1 = tau2*randn(n, nsim); 
  v2 = tau2*randn(n, nsim); 
  
  % Generate Uniform(0,1) random variables in matrix n x nsim x p
  % dimensions, they are used as covariates, p gives us the number of
  % covariates used in the model. Hence the data space is (0,1)^p
  
  x = rand(n, nsim, p);  
  s_2 = sqrt(p); 
  
  % generate the index (linear combinations) with given parameters, this
  % gives us the matrix of dimension n x nsim , we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), now this z can be fed into the Local frechet regression.
  
  z= (b(1)*x(:,:,1) + b(2)*x(:,:,2)+ b(3)*x(:,:,3))./s_2; 
  
  % Map the point in (-1,1) on to the circumference of unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg1 = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
  % defining a set of functions pre-defined on a unit sphere of p-dimensions
  M = spherefactory(3); % define structures fr, ops to store the pre-defined functions on a sphere
  fr = struct(); 
  ops = struct();
  fr.M = M; 
  ops.verbosity = 0;
 
  % create the variables for storing the simulated response data based on
  % our assumptions on theta and covariates: creates a list of nsim number
  % of matrices of size 3 x n, for each simulation we have 
  
  Y_true1 = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  Y_true2 = Y_true1; 
  
  % for each simulation, generate response data by first transporting the single index
  % onto the unit sphere, then map that point to the tangent space, generates a single point by adding noise created by
  % normal samples in u1, u2.
  
  for i=1:nsim
     Y_true1{i} = cell2mat(arrayfun(@(j) add_noise(mreg1{i}(:, j), [u1(j, i) u2(j, i)]), 1:size(z,1), 'UniformOutput', false)); 
  end
  
  for i=1:nsim
     Y_true2{i} = cell2mat(arrayfun(@(j) add_noise(mreg1{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
% Creating the array with necessary dimensions to store MSPE in both high and low noises respectively 
mspe_ln = arrayfun(@(i) zeros(length(h1), size(start,1)), 1:nsim, 'UniformOutput', false);
mspe_hn = arrayfun(@(i) zeros(length(h2), size(start,1)), 1:nsim, 'UniformOutput', false);

pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim % starting the parallel for loop
%for i = 1:nsim
disp(i)
for j=1:size(start,1)
    
    xtemp = [x(:,i,1), x(:,i,2), x(:,i,3)];     % getting the covariates for the i-th simulation

    for l1 = 1:length(h1) % for every bandwidth in the low noise
       disp(l1)
       Y_temp1 = Y_true1{i}; % creating temporary response variable to be input for cost function
       cost_fun = @(eta) nm_cost(Y_temp1, eta, h1(l1), xtemp); % define cost function 
       
       % running the optimization for starting point 
       [eta_opt1, fval] = fminsearchbnd(cost_fun, [start(j,1), start(j,2)], [-pi,-pi], [pi, pi] );
        mspe_ln{i}(l1,j) = fval;            % storing minimum cost function
        eta_opt_mspe_ln{i}(l1,:,j) = eta_opt1;   % storing local minima
       
    end
    
    % repeat the above procedure for high noise scenario
    for l2 = 1:length(h2) % bandwidth range for high noise scenario
     
     disp(l2)
       Y_temp2 = Y_true2{i}; % creating temporary response variable to be input for cost function
       cost_fun = @(eta) nm_cost(Y_temp2, eta, h2(l2), xtemp); % define cost function 
       
       % running the optimization for starting point 
       [eta_opt2, fval] = fminsearchbnd(cost_fun, [start(j,1), start(j,2)], [-pi,-pi],[pi, pi]);
       
       mspe_hn{i}(l2) = fval;          % storing minimum cost function
       eta_opt_mspe_hn{i}(l2,:)= eta_opt2;   % storing local minima
     
    end
  
end

end

delete(pool) % shutting down parallel pool 

% saving the results and output for further analysis 
save Sphere_results_n50_nsim200_p3_neldermead.mat mreg1 * mspe_hn * mspe_ln * eta_opt_mspe_hn * eta_opt_mspe_ln * h1 * h2























%% Sphere simulations n=50, p=4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Adding path to folders
addpath(genpath('manopt'))
addpath(genpath('FMINSEARCHBND'))

% For beta in circumference of a unit circle;

nsim = 200;  % number of simulated datasets
n=50;        % number of observations 
p=4;         % dimension of covariate vector 
s=19878;     % random number generators
rng(s);      % setting seed for reproducible simulations

% Variances for normal data simulation
tau1=0.2;   % Low Noise 
tau2=0.6;  % High Noise

% function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true beta parameter
b =[0.6087, 0.1332, -0.287, -sqrt(1-0.6087^2-0.1332^2-0.287^2)];

% Creating the grid for theta on hypersphere p=3, using polar coordinates
phi = linspace(-pi+(pi/5), pi-(pi/5),3); %Choose the grid for eta to start algorithm  
[S1 S2 S3] = ndgrid(phi', phi', phi');
start = [S1(:), S2(:), S3(:)]; % creating the grid
% writematrix(start,'M.csv') % write the data frame as .csv file for further use 

% plot check
%plot3(start(:,1), start(:,2), start(:,3), 'r.');
%ylim = [-pi, pi];
%xlim = [-pi, pi];
%zlim = [-pi, pi];
%size(start)       % check the size of the 

h1 = exp(log(0.45):0.02:log(0.85));     % bandwidth sequence for low noise data
h2 = exp(log(0.45):0.02:log(0.85));    % bandwidth sequence for high noise data

pm = struct('s', s, 'nsim', nsim, 'n', n, 'h1', h1, 'h2', h2, 'tau1', tau1, 'tau2', tau2);

eta_opt_mspe_hn = arrayfun(@(l) zeros(length(h2),p-1, size(start,1)), 1:nsim, 'UniformOutput', false);
eta_opt_mspe_ln = arrayfun(@(l) zeros(length(h1),p-1, size(start,1)), 1:nsim, 'UniformOutput', false);

% Generating the data from fixed \theta_0 value, here we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

  % Generate N(0, (0.2)^2) random numbers in matrix n x nsim in 2 sets u1 and u2 respectively
  u1 = tau1*randn(n, nsim); 
  u2 = tau1*randn(n, nsim); 
  
  % Generate N(0, (0.35)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  v1 = tau2*randn(n, nsim); 
  v2 = tau2*randn(n, nsim); 
  
  % Generate Uniform(0,1) random variables in matrix n x nsim x p
  % dimensions, they are used as covariates, p gives us the number of
  % covariates used in the model. Hence the data space is (0,1)^p
  
  x = rand(n, nsim, p);  
  s_2 = sqrt(p); 
  
  % generate the index (linear combinations) with given parameters, this
  % gives us the matrix of dimension n x nsim , we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), now this z can be fed into the Local frechet regression.
  
  z = (b(1)*x(:,:,1) + b(2)*x(:,:,2)+ b(3)*x(:,:,3)+ b(4)*x(:,:,4))./s_2; 
  
  % Map the point in (-1,1) on to the circumference of unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg1 = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
  % defining a set of functions pre-defined on a unit sphere of p-dimensions
  M = spherefactory(3); % define structures fr, ops to store the pre-defined functions on a sphere
  fr = struct(); 
  ops = struct();
  fr.M = M; 
  ops.verbosity = 0;
 
  % create the variables for storing the simulated response data based on
  % our assumptions on theta and covariates: creates a list of nsim number
  % of matrices of size 3 x n, for each simulation we have 
  
  Y_true1 = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  Y_true2 = Y_true1; 
  
  % for each simulation, generate response data by first transporting the single index
  % onto the unit sphere, then map that point to the tangent space, generates a single point by adding noise created by
  % normal samples in u1, u2.
  
  for i=1:nsim
     Y_true1{i} = cell2mat(arrayfun(@(j) add_noise(mreg1{i}(:, j), [u1(j, i) u2(j, i)]), 1:size(z,1), 'UniformOutput', false)); 
  end
  
  for i=1:nsim
     Y_true2{i} = cell2mat(arrayfun(@(j) add_noise(mreg1{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
% Creating the array with necessary dimensions to store MSPE in both high and low noises respectively 
mspe_ln = arrayfun(@(i) zeros(length(h1), size(start,1)), 1:nsim, 'UniformOutput', false);
mspe_hn = arrayfun(@(i) zeros(length(h2), size(start,1)), 1:nsim, 'UniformOutput', false);

pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim % starting the parallel for loop
%for i = 1:nsim
disp(i)

for j=1:size(start,1) % for every starting point for algorithm
   
    xtemp = [x(:,i,1), x(:,i,2), x(:,i,3), x(:,i,4)];     % getting the covariates for the i-th simulation

    for l1 = 1:length(h1) % for every bandwidth in the low noise
       disp(l1)
       
       Y_temp1 = Y_true1{i}; % creating temporary response variable to be input for cost function
       cost_fun = @(eta) nm_cost(Y_temp1, eta, h1(l1), xtemp); % define cost function 
       
       % running the optimization for starting point with specified upper
       % lower bounds for all variables.
       [eta_opt1, fval] = fminsearchbnd(cost_fun, [start(j,1),start(j,2),start(j,3)]', [-pi,-pi,-pi],[pi, pi,pi] );
       
       mspe_ln{i}(l,j) = fval;            % storing minimum cost function
       eta_opt_mspe_ln{i}(l,:, j) = eta_opt1;   % storing local minima
       
    end
    
    % we repeat the similar steps above for the high noise scenario~
    
    for l2 = 1:length(h2) % bandwidth range for high noise scenario
       disp(l2)
      
       Y_temp2 = Y_true2{i};  % creating temporary response variable to be input for cost function
       cost_fun = @(eta) nm_cost(Y_temp2, eta, h2(l2), xtemp); % define cost function 
       
       % running the optimization for starting point 
       [eta_opt2, fval] = fminsearchbnd(cost_fun, [start(j,1),start(j,2),start(j,3)]', [-pi,-pi,-pi],[pi, pi,pi]);
       
       mspe_hn{i}(l2,j) = fval;          % storing minimum cost function
       eta_opt_mspe_hn{i}(l2,:,j)= eta_opt2;   % storing local minima
     
    end
  
end

end

delete(pool) % shutting down parallel pool 

% save the output for further analysis
save Sphere_results_n50_nsim200_p4_neldermead.mat mreg1 * start * mspe_hn * mspe_ln * eta_opt_mspe_hn * eta_opt_mspe_ln * h1 * h2

