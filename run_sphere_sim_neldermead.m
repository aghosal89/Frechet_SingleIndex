%% Computation of MSPE/MSEE for various simulation scenarios using Nelder-Mead optimization

%% Sphere simulations n, p=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
% Adding path to 'manopt' folder
addpath(genpath('manopt 2'))
addpath(genpath('FMINSEARCHBND'))

% For beta in circumference of a unit circle;

nsim = 200;   % number of simulations
n=50;        % number of observations 
p=2;         % dimension of theta vector 
s=19878;     % random number generators
rng(s);      % setting seed for random number generator

% Variances for normal data simulation
tau=0.2;    % Noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true theta
b =[cos(-0.704), sin(-0.704)];

% Creating the starting points for eta on hypersphere p=2
nsp = 3; % number of starting points
sc= pi/nsp;
start = linspace(-pi/2+(sc/2),pi/2-sc/2,nsp)';         

% Generating the data from fixed \theta_0 value, here we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 
  
  % Generate N(0, (tau)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  v1 = tau*randn(n, nsim); 
  v2 = tau*randn(n, nsim); 
  
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
  mreg = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
  % defining a set of functions pre-defined on a unit sphere of p-dimensions
  M = spherefactory(3); % define structures fr, ops to store the pre-defined functions on a sphere
  fr = struct(); 
  ops = struct();
  fr.M = M; 
  ops.verbosity = 0;
 
  % create the variables for storing the simulated response data based on
  % our assumptions on theta and covariates: creates a list of nsim number
  % of matrices of size 3 x n, for each simulation we have 
  
  Y_true = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  
  % for each simulation, generate response data by first transporting the single index
  % onto the unit sphere, then map that point to the tangent space, generates a single point by adding noise created by
  % normal samples in v1, v2.
  
  for i=1:nsim
     Y_true{i} = cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
% choosing a range bandwidths~
% Creating the grid for theta on hypersphere p=2, using polar coordinates
nsp = 20;
sc= pi/nsp;
f = linspace(-pi/2+sc/2, pi/2-sc/2, nsp);
beta = [cos(f)', sin(f)'];

bw_min= zeros(size(beta,1),1);
bw_max= zeros(size(beta,1),1);
temp=zeros(nsim,1);
for k=1:size(beta,1)
    for i=1:nsim
        temp = (beta(k,1)*x(:,i,1) + beta(k,2)*x(:,i,2))./sqrt(p);
        tempd=abs(repmat(temp, 1, length(temp)) - repmat(temp', length(temp), 1));
        tempc=sort(tempd,2);
        bw_min(k)= max(tempc(:,4));
        bw_max(k)=max(tempd(:));
    end
end    

%h2=exp(linspace(log(max(bw_min)*1.25),log(max(bw_max)),15))
h =0.437;

pm = struct('s', s, 'nsim', nsim, 'n', n, 'h', h, 'tau', tau);

% Creating the array with necessary dimensions to store MSPE in both high and low noises respectively 

mspe = arrayfun(@(i) zeros(length(h), size(start,1)), 1:nsim, 'UniformOutput', false);

eta_opt = arrayfun(@(l) zeros(length(h), p-1, size(start,1)), 1:nsim, 'UniformOutput', false);

pool=parpool(6);   % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim  % starting the parallel for loop
%for i = 1:nsim
disp(i)

Y_temp = Y_true{i};             % creating temporary response variable to be input for cost function
xtemp = [x(:,i,1), x(:,i,2)];   % getting the covariates for the i-th simulation
              
for j=1:size(start,1)

    for l = 1:length(h) % bandwidth range for high noise scenario
    disp(['NM_nsim200_n50_p2_LN' '_' num2str(i, '%03d') '_' 'start' '_' num2str(j, '%03d') '_' 'bandwidth' '_' num2str(l, '%02d')])
    
    cost_fun = @(eta) nm_cost(Y_temp, eta, h(l), xtemp); % define cost function 
       
    % running the optimization for starting point 
    [eta_opt1, fval] = fminsearchbnd(cost_fun, start(j,:), -pi/2, pi/2);
    mspe{i}(l, j) = fval;           % storing minimum cost function
    eta_opt{i}(l,:,j) = eta_opt1;   % storing local minima
  
    end
  
end

end

delete(pool) % shutting down parallel pool 

save NM_Sphere_results_n50_nsim200_p2_LN.mat * mspe * eta_opt * h


















%% Siulation results n, p=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

% Adding path to files and folder
addpath(genpath('manopt 2'))
addpath(genpath('FMINSEARCHBND'))

% For beta in circumference of a unit circle;

nsim = 200;   % number of simulations
n=200;        % number of observations 
p=3;         % dimension of theta vector 
s=19878;     % random number generators
rng(s);      % setting seed for random number generator

%Variances for normal data simulation
tau=0.6;  % Noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true beta parameter
b=[0.386, 0.181, -sqrt(1-0.386^2-0.181^2)];

% Creating the grid for theta on hypersphere p=3, using polar coordinates
nsp = 4;
sc=pi/nsp;
phi = linspace(-pi/2+(sc/2), pi/2-(sc/2),nsp); %Choose the grid for eta 
[S1 S2] = ndgrid(phi', phi');
start = [S1(:), S2(:)]; % creating the grid

% Generating the data from fixed \theta_0 value, here we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

  % Generate N(0, (tau)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  v1 = tau*randn(n, nsim); 
  v2 = tau*randn(n, nsim); 
  
  % Generate Uniform(0,1) random variables in matrix n x nsim x p
  % dimensions, they are used as covariates, p gives us the number of
  % covariates used in the model. Hence the data space is (0,1)^p
  
  x = rand(n, nsim, p);  
  
  % generate the index (linear combinations) with given parameters, this
  % gives us the matrix of dimension n x nsim , we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), now this z can be fed into the Local frechet regression.
  
  z= (b(1)*x(:,:,1) + b(2)*x(:,:,2)+ b(3)*x(:,:,3))./sqrt(p); 
  
  % Map the point in (-1,1) on to the circumference of unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
  % defining a set of functions pre-defined on a unit sphere of p-dimensions
  M = spherefactory(3); % define structures fr, ops to store the pre-defined functions on a sphere
  fr = struct(); 
  ops = struct();
  fr.M = M; 
  ops.verbosity = 0;
 
  % create the variables for storing the simulated response data based on
  % our assumptions on theta and covariates: creates a list of nsim number
  % of matrices of size 3 x n, for each simulation we have 
  
  Y_true = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  
  % for each simulation, generate response data by first transporting the single index
  % onto the unit sphere, then map that point to the tangent space, generates a single point by adding noise created by
  % normal samples in v1, v2.
  
  for i=1:nsim
     Y_true{i} = cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
% choosing a range bandwidths~
% Creating the grid for theta on hypersphere p=2, using polar coordinates
nsp=30;
sc= pi/nsp;
f=linspace(-pi/2+sc/2,pi/2-sc/2, 30);
[S1 S2] = ndgrid(f', f');
X1= cos(S1).*cos(S2);
X2= cos(S1).*sin(S2);
X3= sin(S1);

beta = [X1(:), X2(:), X3(:)];
beta = beta(beta(:,1)>0.01, :);

bw_min= zeros(size(beta,1),1);
bw_max= zeros(size(beta,1),1);
temp=zeros(nsim,1);
for k=1:size(beta,1)
    for i=1:nsim
        temp = (beta(k,1)*x(:,i,1) + beta(k,2)*x(:,i,2)+ beta(k,3)*x(:,i,3))./sqrt(p);
        tempd=abs(repmat(temp, 1, length(temp)) - repmat(temp', length(temp), 1));
        tempc=sort(tempd,2);
        bw_min(k)= max(tempc(:,4));
        bw_max(k)= max(tempd(:));
    end
end    

%h=exp(linspace(log(max(bw_min)*1.25),log(max(bw_max)),10));
h= .374; 

pm = struct('s', s, 'nsim', nsim, 'n', n, 'h', h, 'tau', tau);

% Creating the array with necessary dimensions to store MSPE in both high and low noises respectively 

mspe = arrayfun(@(i) zeros(length(h), size(start,1)), 1:nsim, 'UniformOutput', false);

eta_opt = arrayfun(@(l) zeros(length(h), p-1, size(start,1)), 1:nsim, 'UniformOutput', false);

pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim  % starting the parallel for loop
%for i = 1:nsim
disp(i)
for j=1:size(start,1)
    
    xtemp = [x(:,i,1), x(:,i,2), x(:,i,3)];     % getting the covariates for the i-th simulation

    % estimation procedure for high noise scenario
    for l2 = 1:length(h) % bandwidth range for high noise scenario
       disp(['NM_nsim200_n200_p3_HN' '_' num2str(i, '%03d') '_' 'start' '_' num2str(j, '%03d') '_' 'bandwidth' '_' num2str(l2, '%02d')])

       Y_temp = Y_true{i}; % creating temporary response variable to be input for cost function
       cost_fun = @(eta) nm_cost(Y_temp, eta, h(l2), xtemp); % define cost function 
       
       % running the optimization for starting point 
       [eta_opt2, fval] = fminsearchbnd(cost_fun, start(j,:), [-pi/2,-pi/2],[pi/2, pi/2]);
       
       mspe{i}(l2,j) = fval;          % storing minimum cost function
       eta_opt{i}(l2,:,j)= eta_opt2;   % storing local minima
     
    end
  
end

end

delete(pool) % shutting down parallel pool 

% saving the results and output for further analysis 
save NM_Sphere_results_n200_nsim200_p3_HN.mat * mspe * eta_opt * h















%% Sphere simulations n, p=4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%Adding path to files and folder
addpath(genpath('manopt 2'))
addpath(genpath('FMINSEARCHBND'))

% For beta in circumference of a unit circle;

nsim = 200;  % number of simulated datasets
n=100;        % number of observations 
p=4;         % dimension of covariate vector 
s=19878;     % random number generators
rng(s);      % setting seed for reproducible simulations

% Variances for normal data simulation
tau=0.2;  % Noise

% function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true beta parameter
b =[0.6087, 0.1332, -0.287, -sqrt(1-0.6087^2-0.1332^2-0.287^2)];

% Creating the grid for theta on hypersphere p=3, using polar coordinates
nsp = 3; % number of starting points for each dimention of polar coordinates
sc= pi/nsp;
phi = linspace(-pi/2+(sc/2), pi/2-(sc/2),nsp); %Choose the grid for eta to start algorithm  
[S1 S2 S3] = ndgrid(phi', phi', phi');
start = [S1(:), S2(:), S3(:)]; % creating the grid

% Generating the data from fixed \theta_0 value, here we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

  % Generate N(0, (tau)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  v1 = tau*randn(n, nsim); 
  v2 = tau*randn(n, nsim); 
  
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
  mreg = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
  % defining a set of functions pre-defined on a unit sphere of p-dimensions
  M = spherefactory(3); % define structures fr, ops to store the pre-defined functions on a sphere
  fr = struct(); 
  ops = struct();
  fr.M = M; 
  ops.verbosity = 0;
 
  % create the variables for storing the simulated response data based on
  % our assumptions on theta and covariates: creates a list of nsim number
  % of matrices of size 3 x n, for each simulation we have 
  
  Y_true = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  
  % for each simulation, generate response data by first transporting the single index
  % onto the unit sphere, then map that point to the tangent space, generates a single point by adding noise created by
  % normal samples in u1, u2.
  
  for i=1:nsim
     Y_true{i} = cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
% choosing a range bandwidths~
nsp = 20;
sc= pi/nsp;
f=linspace(-pi/2+sc/2,pi/2-sc/2, nsp);
[S1 S2 S3] = ndgrid(f', f', f');
X1= cos(S1).*cos(S2).*cos(S3);
X2= cos(S1).*cos(S2).*sin(S3);
X3= cos(S1).*sin(S2);
X4= sin(S1);

beta = [X1(:), X2(:), X3(:), X4(:)];
beta = beta(beta(:,1)>0.01, :);

bw_min= zeros(size(beta,1),1);
bw_max= zeros(size(beta,1),1);
temp=zeros(nsim,1);
for k=1:size(beta,1)
    for i=1:nsim
        temp = (beta(k,1)*x(:,i,1) + beta(k,2)*x(:,i,2)+ beta(k,3)*x(:,i,3)+beta(k,4)*x(:,i,4))./sqrt(p);
        tempd=abs(repmat(temp, 1, length(temp)) - repmat(temp', length(temp), 1));
        tempc=sort(tempd,2);
        bw_min(k)= max(tempc(:,4));
        bw_max(k)= max(tempd(:));
    end
end    

h=exp(linspace(log(max(bw_min)*1.25),log(max(bw_max)),10));

pm = struct('s', s, 'nsim', nsim, 'n', n, 'h', h, 'tau', tau);

% Creating the array with necessary dimensions to store MSPE in both high and low noises respectively 

mspe = arrayfun(@(i) zeros(length(h), size(start,1)), 1:nsim, 'UniformOutput', false);

eta_opt = arrayfun(@(l) zeros(length(h), p-1, size(start,1)), 1:nsim, 'UniformOutput', false);
  
pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim % starting the parallel for loop
%for i = 1:nsim
disp(i)

for j=1:size(start,1) % for every starting point for algorithm
   
    xtemp = [x(:,i,1), x(:,i,2), x(:,i,3), x(:,i,4)];     % getting the covariates for the i-th simulation

    % we repeat the similar steps above for the high noise scenario~
    
    for l2 = 1:length(h) % bandwidth range for high noise scenario
       disp(['NM_n50_nsim50_p4_LN_nsim' '_' num2str(i, '%03d') '_' 'start' '_' num2str(j, '%03d') '_' 'bandwidth' '_' num2str(l2, '%02d')])
    
       Y_temp = Y_true{i};  % creating temporary response variable to be input for cost function
       cost_fun = @(eta) nm_cost(Y_temp, eta, h(l2), xtemp); % define cost function 
       
       % running the optimization for starting point 
       [eta_opt2, fval] = fminsearchbnd(cost_fun, start(j,:), [-pi/2,-pi/2,-pi/2],[pi/2, pi/2,pi/2]);
       
       mspe{i}(l2,j) = fval;           % storing minimum cost function
       eta_opt{i}(l2,:,j)= eta_opt2;   % storing local minima
     
    end
  
end

end

delete(pool) % shutting down parallel pool 

% save the output for further analysis
save NM_Sphere_results_n50_nsim200_p4_HN.mat * mspe * eta_opt * h




