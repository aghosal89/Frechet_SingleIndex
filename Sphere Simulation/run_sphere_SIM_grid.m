%% For computing Mean Square Prediction Error (MSPE), Mean Square Estimation Error (MSEE)

%% Sphere simulations n=50, p=4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;   % clearning all working environment
% generating path for contents of folder 'Manopt_6.0'
addpath(genpath('manopt'))

nsim = 200;   % number of simulations
n=200;        % number of observations 
p=4;         % dimension of theta vector 
s=19378;     % random number generators
rng(s);      % setting seed for random number generator

%Variances for normal data simulation
tau=0.6;  % High Noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true theta for single index model~
b =[0.6087, 0.1332, -0.287, -sqrt(1-0.6087^2-0.1332^2-0.287^2)];

% Creating the grid for theta on hypersphere p=4, using polar coordinates
nsp=10;
spc= pi/nsp;
f = linspace(-pi/2+spc/2,pi/2-spc/2,((pi-spc)/spc)+1); %choosing the increments
[S1 S2 S3] = ndgrid(f', f', f'); % creating the grid
X = [S1(:), S2(:), S3(:)];

beta=zeros(size(X,1), size(X,2)+1);
for i=1:size(beta,1)
    beta(i,:)= polar2cart(X(i,:),1);
end    

% Generating the data from fixed theta value, here we assume that the
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
  % gives us the matrix of dimension n x nsim, we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), now this z can be fed into the Local frechet regression.
  
  z= (b(1)*x(:,:,1) + b(2)*x(:,:,2) + b(3)*x(:,:,3)+ b(4)*x(:,:,4))./sqrt(p); 
  
  % Map the point in (-1,1) on to the circumference of unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
  % Grids and variable storage
  
  %defining a set of functions pre-defined on a unit sphere of p-dimensions
  M = spherefactory(p); % define structures fr, ops to store the pre-defined functions on a sphere
  fr = struct(); 
  ops = struct();
  fr.M = M; 
  ops.verbosity = 0;
 
  % create the variables for storing the simulated response data based on
  % our assumptions on theta and covariates: creates a list of nsim number
  % of matrices of size 3 x n, for each simulation we have 
  
  Y = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  
  % for each simulation, generate response data by first transporting the single index
  % onto the unit sphere, then map that point to the tangent space, generates a single point by 
  % adding noise created by normal samples in v1, v2.
   
  for i=1:nsim
     Y{i} = cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end

% Creating the grid for theta on hypersphere, using polar coordinates
bw_min= zeros(size(beta,1),nsim);
bw_max= zeros(size(beta,1),nsim);

for k=1:size(beta,1)
    for i=1:nsim
        temp = (beta(k,1)*x(:,i,1) + beta(k,2)*x(:,i,2)+ beta(k,3)*x(:,i,3)+ beta(k,4)*x(:,i,4))./sqrt(p);
        tempd=abs(repmat(temp, 1, length(temp)) - repmat(temp', length(temp), 1));
        tempc=sort(tempd,2);
        bw_min(k,i)= max(tempc(:,4));
        bw_max(k,i)=max(tempd(:));
    end
end    

h=exp(linspace(log(max(bw_min(:))*1.25),log(max(bw_max(:))),10));

% Creating the array with necessary dimensions to store MSPE and MSEE respectively
mspe = arrayfun(@(i) zeros(length(beta), length(h)), 1:nsim, 'UniformOutput', false);
msee = arrayfun(@(i) zeros(length(beta), length(h)), 1:nsim, 'UniformOutput', false);

pm = struct('s', s, 'nsim', nsim, 'n', n, 'tau', tau, 'beta', beta);

pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim % starting the parallel for loop
% for i = 1:nsim
disp(i)

for r = 1:size(beta,1) % starting for loop on the theta parameter.
    
disp(r)
xtemp = [x(:,i,1), x(:,i,2), x(:,i,3), x(:,i,4)];     % getting the covariates for the i-th simulation
ztemp = (xtemp*beta(r,:)')./sqrt(p);  % dividing by sqrt(p) so that the entries in ztemp in [-1,1].
  
    for l = 1:length(h)

       Y_fr = get_sphere_fit_FSI(Y{i}, ztemp, h(l)); % getting the frechet mean for given theta, covariates
       mspe{i}(r,l) = mean((arrayfun(@(g) acos(Y{i}(:, g)'*Y_fr(:,g)),1:n)).^2); % compute the MSPE
       msee{i}(r,l) = mean((arrayfun(@(g) acos(mreg{i}(:, g)'*Y_fr(:,g)),1:n)).^2);    % compute the MSEE

    end
  
end

end
delete(pool) % shutting down parallel pool 

mspe_mat = NaN(size(beta,1), nsim, length(h));

for i=1:nsim
    for r=1:length(beta)
        for q=1:length(h)
            mspe_mat(r, i, q) = (mspe{i}(r,q));
        end
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_sum_mat = mspe{1};

for i=2:nsim
    mspe_sum_mat = mspe_sum_mat + mspe{i}; 
end

mspe_ave_mat = mspe_sum_mat./nsim;

[~, c2] = find(mspe_ave_mat==min(mspe_ave_mat(:)));

% estimate of bandwidth
h_opt= h(c2);

beta_opt =zeros(nsim, size(beta,2));
for i=1:nsim
    beta_opt(i,:)= beta(find(mspe_mat(:,i,c2)==min(mspe_mat(:,i,c2))),:);
end    

% mean and standard error of the theta estimates
mean(acos(beta_opt*b')) 
std(acos(beta_opt*b'))

save Grid_Sphere_results_n200_nsim200_p4_HN.mat mreg * mspe * msee * h * Y















%% Sphere simulations n, p=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clearing the working environment
clear all;
% adding path to the folder
addpath(genpath('manopt'))

% For beta in circumference of a unit circle;

nsim = 200;  % number of simulations
n=100;        % number of observations 
p=3;         % dimension of theta vector 
s=19379;     % random number generators
rng(s);      % setting seed for random number generator

%Variances for normal data simulation
%tau=0.6;  % High Noise
tau=0.2;  % Low Noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true beta parameter
b=[0.386, 0.181, -sqrt(1-0.386^2-0.181^2)];

% Generating the data from fixed theta value, here we assume that the
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
  % gives us the matrix of dimension n x nsim, we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), now this z can be fed into the Local frechet regression.
  
  z= (b(1)*x(:,:,1) + b(2)*x(:,:,2) + b(3)*x(:,:,3))./sqrt(p); 
  
  % Map the point in (-1,1) on to the circumference of unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
  % Grids and variable storage
  
  %defining a set of functions pre-defined on a unit sphere of p-dimensions
  M = spherefactory(p); % define structures fr, ops to store the pre-defined functions on a sphere
  fr = struct(); 
  ops = struct();
  fr.M = M; 
  ops.verbosity = 0;
 
  % create the variables for storing the simulated response data based on
  % our assumptions on theta and covariates: creates a list of nsim number
  % of matrices of size 3 x n, for each simulation we have 
  
  Y = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  
  % for each simulation, generate response data by first transporting the single index
  % onto the unit sphere, then map that point to the tangent space, generates a single point by adding noise created by
  % normal samples in v1, v2.
    
  for i=1:nsim
     Y{i} = cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
% Creating the grid for theta on hypersphere p=3, using polar coordinates
nsp=20;
spc= pi/nsp;
f = linspace(-pi/2+spc/2,pi/2-spc/2,((pi-spc)/spc)+1); %choosing the increments
[S1 S2] = ndgrid(f', f'); % creating the grid
X = [S1(:), S2(:)];

beta=zeros(size(X,1), size(X,2)+1);
for i=1:size(beta,1)
    beta(i,:)= polar2cart(X(i,:),1);
end    

% creating a grid of bandwidths
bw_min= zeros(size(beta,1),nsim);
bw_max= zeros(size(beta,1),nsim);

for k=1:size(beta,1)
    for i=1:nsim
        temp = (beta(k,1)*x(:,i,1) + beta(k,2)*x(:,i,2)+ beta(k,3)*x(:,i,3))./sqrt(p);
        tempd=abs(repmat(temp, 1, length(temp)) - repmat(temp', length(temp), 1));
        tempc=sort(tempd,2);
        bw_min(k,i)= max(tempc(:,4));
        bw_max(k,i)=max(tempd(:));
    end
end    

h=exp(linspace(log(max(bw_min(:))*1.25),log(max(bw_max(:))),10));

% Creating the array with necessary dimensions to store MSPE, MSEE in high noise respectively 
mspe = arrayfun(@(i) zeros(length(beta), length(h)), 1:nsim, 'UniformOutput', false);
msee = arrayfun(@(i) zeros(length(beta), length(h)), 1:nsim, 'UniformOutput', false);

pm = struct('s', s, 'nsim', nsim, 'n', n, 'tau', tau, 'beta', beta);

pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim % starting the parallel for loop
%for i = 1:nsim

for r = 1:length(beta) % starting for loop on the theta parameter.

xtemp = [x(:,i,1), x(:,i,2), x(:,i,3)];     % getting the covariates for the i-th simulation
ztemp = (xtemp*beta(r,:)')./sqrt(p);  % dividing by sqrt(p) so that the entries in ztemp in [-1,1].

    for l = 1:length(h)
       disp(['grid_n50nsim50p3LN_nsim' '_' num2str(i, '%03d') '_' 'beta' '_' num2str(r, '%03d') '_' 'bandwidth' '_' num2str(l, '%02d')])
       Y_fr = get_sphere_fit_FSI(Y{i}, ztemp, h(l)); % getting the frechet mean for given theta, covariates
       mspe{i}(r,l) = mean((arrayfun(@(g) acos(Y{i}(:, g)'*Y_fr(:,g)),1:n)).^2); % compute the MSPE
       msee{i}(r,l) = mean((arrayfun(@(g) acos(mreg{i}(:, g)'*Y_fr(:,g)),1:n)).^2);    % compute the MSEE

    end
  
end

end
delete(pool) % shutting down parallel pool 

mspe_mat = NaN(size(beta,1), nsim, length(h));

for i=1:nsim
    for r=1:length(beta)
        for q=1:length(h)
            mspe_mat(r, i, q) = (mspe{i}(r,q));
        end
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_sum_mat = mspe{1};

for i=2:nsim
    mspe_sum_mat = mspe_sum_mat + mspe{i}; 
end

mspe_ave_mat = mspe_sum_mat./nsim;

[~, c2] = find(mspe_ave_mat==min(mspe_ave_mat(:)));

% estimate of bandwidth
h_opt= h(c2);

beta_opt =zeros(nsim, size(beta,2));
for i=1:nsim
    beta_opt(i,:)= beta(find(mspe_mat(:,i,c2)==min(mspe_mat(:,i,c2))),:);
end    

% mean and standard error of the theta estimates
mean(acos(beta_opt*b')) 
std(acos(beta_opt*b'))

save Grid_Sphere_results_n50_nsim50_p3_HN.mat * 
















%% Sphere simulations n, p=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; % clearning working environment~
%Adding path to files and folder
addpath(genpath('manopt'))

nsim = 200;  % number of simulations
n=50;        % number of observations 
p=2;         % dimension of theta vector 
s=19380;     % random number generators
rng(s);      % setting seed for random number generator

%Variances for normal data simulation
tau =0.2; % noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true theta
b =[cos(-0.704), sin(-0.704)];

% Creating the grid for theta on hypersphere p=4, using polar coordinates
nsp=100;
spc= pi/nsp;
f = linspace(-pi/2+spc/2,pi/2-spc/2,((pi-spc)/spc)+1); %choosing the increments
[S1] = ndgrid(f'); % creating the grid
X = [S1(:)];

beta=zeros(size(X,1), size(X,2)+1);
for i=1:size(beta,1)
    beta(i,:)= polar2cart(X(i,:),1);
end    

% Generating the data from fixed theta value, here we assume that the
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
 % gives us the matrix of dimension n x nsim, we divide by sqrt(p) so that all elements of z 
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
  
  Y = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  
  % for each simulation, generate response data by first transporting the single index
  % onto the unit sphere, then map that point to the tangent space, generates a single point by adding noise created by
  % normal samples in v1, v2.
  
  for i=1:nsim
     Y{i} = cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
% Creating the grid for theta on hypersphere p=3, using polar coordinates
bw_min= zeros(size(beta,1),nsim);
bw_max= zeros(size(beta,1),nsim);

for k=1:size(beta,1)
    for i=1:nsim
        temp = (beta(k,1)*x(:,i,1) + beta(k,2)*x(:,i,2))./sqrt(p);
        tempd=abs(repmat(temp, 1, length(temp)) - repmat(temp', length(temp), 1));
        tempc=sort(tempd,2);
        bw_min(k,i)= max(tempc(:,4));
        bw_max(k,i)=max(tempd(:));
    end
end    

h = exp(linspace(log(max(bw_min(:))*1.25),log(max(bw_max(:))),10));

pm = struct('s', s, 'nsim', nsim, 'n', n, 'tau', tau, 'beta', beta, 'h', h);

% Creating the array with necessary dimensions to store MSPE in both high and low noises respectively 
mspe = arrayfun(@(i) zeros(length(beta), length(h)), 1:nsim, 'UniformOutput', false);
msee = arrayfun(@(i) zeros(length(beta), length(h)), 1:nsim, 'UniformOutput', false);

pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim % starting the parallel for loop
%for i = 1:nsim
xtemp = [x(:,i,1), x(:,i,2)];     % getting the covariates for the i-th simulation

for r = 1:length(beta) % starting for loop on the theta parameter.
    
ztemp = (xtemp*beta(r,:)')./sqrt(p);  % dividing by sqrt(p) so that the entries in ztemp in [-1,1].

    for l = 1:length(h)
       disp(['grid_nsim_n50nsim200p2LN' '_' num2str(i, '%03d') '_' 'beta' '_' num2str(r, '%03d') '_' 'bandwidth' '_' num2str(l, '%02d')])
       
       Y_fr = get_sphere_fit_FSI(Y{i}, ztemp, h(l)); % getting the frechet mean for given theta, covariates
       mspe{i}(r,l) = mean((arrayfun(@(g) acos(Y{i}(:, g)'*Y_fr(:,g)),1:n)).^2); % compute the MSPE
       msee{i}(r,l) = mean((arrayfun(@(g) acos(mreg{i}(:, g)'*Y_fr(:,g)),1:n)).^2);    % compute the MSEE

    end
  
end

end
delete(pool) % shutting down parallel pool 

mspe_mat = NaN(size(beta,1), nsim, length(h));

for i=1:nsim
    for r=1:length(beta)
        for q=1:length(h)
            mspe_mat(r, i, q) = (mspe{i}(r,q));
        end
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_sum_mat = mspe{1};

for i=2:nsim
    mspe_sum_mat = mspe_sum_mat + mspe{i}; 
end

mspe_ave_mat = mspe_sum_mat./nsim;

[~, c2] = find(mspe_ave_mat==min(mspe_ave_mat(:)));

% estimate of bandwidth
h_opt= h(c2);

beta_opt =zeros(nsim, size(beta,2));
for i=1:nsim
    beta_opt(i,:)= beta(find(mspe_mat(:,i,c2)==min(mspe_mat(:,i,c2))),:);
end    

% mean and standard error of the theta estimates
mean(acos(beta_opt*b')) 
std(acos(beta_opt*b'))

save Grid_Sphere_NM_results_n50_nsim200_p2_LN.mat * 


