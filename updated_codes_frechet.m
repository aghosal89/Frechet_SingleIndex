%% For computing Mean Square Prediction Error (MSPE), Mean Square Estimation Error (MSEE)

%% Sphere simulations n=50, p=4, High Noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generating path for files and folders
addpath(genpath('manopt'))

nsim =2;   % number of simulations
n=50;        % number of observations 
p=4;         % dimension of theta vector 
s=19878;     % random number generators
rng(s);      % setting seed for random number generator

%Variances for normal data simulation
tau1=0.2;   % Low noise 
tau2=0.6;  % High Noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true theta for single index model~
b =[0.6087, 0.1332, -0.287, -sqrt(1-0.6087^2-0.1332^2-0.287^2)];

% Creating the grid for theta on hypersphere p=4, using polar coordinates
phi = 0:(pi/6):(2*pi); %choosing the increments
[S1 S2 S3] = ndgrid(phi', phi', phi'); % creating the grid
X = [S1(:), S2(:), S3(:)];

beta=zeros(size(X,1), size(X,2)+1);
for i=1:size(beta,1)
    beta(i,:)= polar2cart(X(i,:),1);
end    

beta = beta(beta(:,1)>0.02,:); % choosing the first coordinate to be positive
beta = beta(beta(:,1)<0.97,:); % choosing the first coordinate off the boundary

% Creating the grid for theta on hypersphere p=3, using polar coordinates

h1=exp(log(0.4):0.02:log(0.75));   % possible bandwidths on a grid for low noise data
h2=exp(log(0.4):0.02:log(1));

pm = struct('s', s, 'nsim', nsim, 'n', n, 'tau1', tau1,'tau2', tau2,'h2', h2, 'beta', beta);

% Generating the data from fixed theta value, here we assume that the
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
  % gives us the matrix of dimension n x nsim, we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), now this z can be fed into the Local frechet regression.
  
  z= (b(1)*x(:,:,1) + b(2)*x(:,:,2) + b(3)*x(:,:,3)+ b(4)*x(:,:,4))./sqrt(p); 
  
  % Map the point in (-1,1) on to the circumference of unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg1 = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
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
mspe_ln = arrayfun(@(i) zeros(length(beta), length(h1)), 1:nsim, 'UniformOutput', false);
mspe_hn = arrayfun(@(i) zeros(length(beta), length(h2)), 1:nsim, 'UniformOutput', false);

% Creating the array with necessary dimensions to store MSEE in both high and low noises respectively
msee_ln = arrayfun(@(i) zeros(length(beta), length(h1)), 1:nsim, 'UniformOutput', false);
msee_hn = arrayfun(@(i) zeros(length(beta), length(h2)), 1:nsim, 'UniformOutput', false);


pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim % starting the parallel for loop
% for i = 1:nsim
disp(i)

for r = 1:size(beta,1) % starting for loop on the theta parameter.
    
disp(r)
xtemp = [x(:,i,1), x(:,i,2), x(:,i,3), x(:,i,4)];     % getting the covariates for the i-th simulation
ztemp = (xtemp*beta(r,:)')./sqrt(p);  % dividing by sqrt(p) so that the entries in ztemp in [-1,1].

    for l = 1:length(h1)
       Y_fr_lo = get_sphere_fit_SIM2(Y_true1{i}, ztemp, h1(l)); % getting the frechet mean for given theta, covariates
       mspe_ln{i}(r,l) = mean((arrayfun(@(g) acos(Y_true1{i}(:, g)'*Y_fr_lo(:,g)),1:n)).^2); % compute the MSPE
       msee_ln{i}(r,l) = mean((arrayfun(@(g) acos(mreg1{i}(:, g)'*Y_fr_lo(:,g)),1:n)).^2);  % compute the MSEE
    end
    
    for l2 = 1:length(h2)

       Y_fr_hi = get_sphere_fit_SIM2(Y_true2{i}, ztemp, h2(l2)); % getting the frechet mean for given theta, covariates
       mspe_hn{i}(r,l2) = mean((arrayfun(@(g) acos(Y_true2{i}(:, g)'*Y_fr_hi(:,g)),1:n)).^2); % compute the MSPE
       msee_hn{i}(r,l2) = mean((arrayfun(@(g) acos(mreg1{i}(:, g)'*Y_fr_hi(:,g)),1:n)).^2);    % compute the MSEE

    end
  
end

end
delete(pool) % shutting down parallel pool 

save Sphere_results_n50_nsim2_p4.mat mreg1 * mspe_hn * msee_hn * mspe_ln * msee_ln * h1 * h2 * Y_true2 * Y_true1

%the following are some computation necessary for further analyses~

mspe_mat_ln = zeros(size(beta,1), nsim, length(h1));
msee_mat_ln = zeros(size(beta,1), nsim, length(h1));
mspe_mat_hn = zeros(size(beta,1), nsim, length(h2));
msee_mat_hn = zeros(size(beta,1), nsim, length(h2));

% initializing sum over nsim simulations~
mspe_matsum_ln = real(mspe_ln{1});
msee_matsum_ln = real(msee_ln{1});
mspe_matsum_hn = real(mspe_hn{1});
msee_matsum_hn = real(msee_hn{1});

for i=1:nsim
    for r=1:length(beta)
        for l=1:length(h1)
            mspe_mat_ln(r, i, l) = real(mspe_ln{i}(r,l));
            msee_mat_ln(r, i, l) = real(msee_ln{i}(r,l));
        end
        for l=1:length(h2)
            mspe_mat_hn(r, i, l) = real(mspe_hn{i}(r,l));
            msee_mat_hn(r, i, l) = real(msee_hn{i}(r,l));
        end
    end
end

% finding out sum and average of MSPE, MSEE over simulations for given theta, bandwidth  
for i=2:nsim
    
    mspe_matsum_hn = mspe_matsum_hn + real(mspe_hn{i});
    msee_matsum_hn = msee_matsum_hn + real(msee_hn{i});
    mspe_matsum_ln = mspe_matsum_ln + real(mspe_ln{i});
    msee_matsum_ln = msee_matsum_ln + real(msee_ln{i});
    
end

% averaging the MSPE/MSEE values over all nsim simulations
mspe_matsum_hn = mspe_matsum_hn./nsim;
msee_matsum_hn = msee_matsum_hn./nsim;
mspe_matsum_ln = mspe_matsum_ln./nsim;
msee_matsum_ln = msee_matsum_ln./nsim;

[x,y]=find(mspe_matsum_hn == min(mspe_matsum_hn(:)));
mspe_matsum_hn(x,y)
betaopt_mspe_hn= beta(x,:); % optimum theta that minimizes MSPE
h_opt_mspe_hn = h2(y); % optimum bandwidth that minimizes MSPE

[x,y]=find(msee_matsum_hn == min(msee_matsum_hn(:)));
msee_matsum_hn(x,y)
betaopt_msee_hn=  beta(x,:); % optimum theta that minimizes MSEE
h_opt_msee_hn= h2(y); % optimum bandwidth that minimizes MSEE

[x,y]=find(mspe_matsum_ln == min(mspe_matsum_ln(:)));
mspe_matsum_ln(x,y)
betaopt_mspe_ln= beta(x,:); % optimum theta that minimizes MSPE
h_opt_mspe_ln= h1(y); % optimum bandwidth that minimizes MSPE

[x,y]=find(msee_matsum_ln == min(msee_matsum_ln(:)));
msee_matsum_ln(x,y)
betaopt_msee_ln = beta(x,:); % optimum theta that minimizes MSEE
h_opt_msee_ln = h1(y); % optimum bandwidth that minimizes MSEE

distv= zeros(nsim, p);
rw1 = find(h1==h_opt_msee_ln);
rw2 = find(h2==h_opt_msee_hn);

for i=1:nsim
    distv(i, :)= beta(find(msee_mat_hn(:,i,rw2)==min(msee_mat_hn(:,i,rw2))),:);
end

distvs = distv(:,1:2);


















%% Sphere simulations n=50, p=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adding path to the folders
addpath(genpath('manopt'))

% For beta in circumference of a unit circle;

nsim = 188;  % number of simulations
n=50;        % number of observations 
p=3;         % dimension of theta vector 
s=19878;     % random number generators
rng(s);      % setting seed for random number generator

%Variances for normal data simulation
tau1=0.2;   % Low noise 
tau2=0.6;  % High Noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true beta parameter
b=[0.386, 0.181, -sqrt(1-0.386^2-0.181^2)];

% Creating the grid for theta on hypersphere p=3, using polar coordinates
phi = 0:(pi/6):(2*pi); %Choose the grid for eta in our document

[S1 S2] = ndgrid(phi', phi');
X1 = cos(S2).*cos(S1); % creating first coordinate
X2 = cos(S2).*sin(S1); % creating second coordinate
X3 = sin(S2);          % creating third coordinate

beta = [X1(:), X2(:), X3(:)]; % creating the grid

beta = beta(beta(:,1)>0.02,:);   %choosing the first coordinate to be positive
beta = beta(beta(:,1)<0.98,:);   %choosing the first coordinate not on boundary

size(beta)

h1= exp(log(0.46):0.02:log(0.7));   % possible bandwidths on a grid for low noise data
h2= exp(log(0.5):0.02:log(0.75));   % possible bandwidths on a grid for high noise data
pm = struct('s', s, 'nsim', nsim, 'n', n, 'tau1', tau1,'tau2', tau2, 'h2', h2, 'beta', beta);

% Generating the data from fixed theta value, here we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

  % Generate N(0, (0.2)^2) random numbers in matrix n x nsim in 2 sets u1 and u2 respectively
  u1 = tau1*randn(n, nsim); 
  u2 = tau1*randn(n, nsim); 
  
  % Generate N(0, (0.6)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  v1 = tau2*randn(n, nsim); 
  v2 = tau2*randn(n, nsim); 
  
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
  mreg1 = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
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
mspe_ln = arrayfun(@(i) zeros(length(beta), length(h1)), 1:nsim, 'UniformOutput', false);
mspe_hn = arrayfun(@(i) zeros(length(beta), length(h2)), 1:nsim, 'UniformOutput', false);

% Creating the array with necessary dimensions to store MSEE in both high and low noises respectively
msee_ln = arrayfun(@(i) zeros(length(beta), length(h1)), 1:nsim, 'UniformOutput', false);
msee_hn = arrayfun(@(i) zeros(length(beta), length(h2)), 1:nsim, 'UniformOutput', false);

pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim % starting the parallel for loop
%for i = 1:nsim
disp(i)

for r = 1:length(beta) % starting for loop on the theta parameter.
    
disp(r)
xtemp = [x(:,i,1), x(:,i,2), x(:,i,3)];     % getting the covariates for the i-th simulation
ztemp = (xtemp*beta(r,:)')./sqrt(p);  % dividing by sqrt(p) so that the entries in ztemp in [-1,1].

    for l = 1:length(h1)

       Y_fr_lo = get_sphere_fit_SIM2(Y_true1{i}, ztemp, h1(l)); % getting the frechet mean for given theta, covariates
       mspe_ln{i}(r,l) = mean((arrayfun(@(g) acos(Y_true1{i}(:, g)'*Y_fr_lo(:,g)),1:n)).^2); % compute the MSPE
       msee_ln{i}(r,l) = mean((arrayfun(@(g) acos(mreg1{i}(:, g)'*Y_fr_lo(:,g)),1:n)).^2);  % compute the MSEE

    end
    
    for l2 = 1:length(h2)

       Y_fr_hi = get_sphere_fit_SIM2(Y_true2{i}, ztemp, h2(l2)); % getting the frechet mean for given theta, covariates
       mspe_hn{i}(r,l2) = mean((arrayfun(@(g) acos(Y_true2{i}(:, g)'*Y_fr_hi(:,g)),1:n)).^2); % compute the MSPE
       msee_hn{i}(r,l2) = mean((arrayfun(@(g) acos(mreg1{i}(:, g)'*Y_fr_hi(:,g)),1:n)).^2);    % compute the MSEE

    end
  
end

end
delete(pool) % shutting down parallel pool 

save Sphere_results_n50_nsim2_p3.mat mreg1 * h2 * mspe_hn * msee_hn * Y_true2

%the following are some computation necessary for further analyses~

mspe_mat_ln = zeros(size(beta,1), nsim, length(h1));
msee_mat_ln = zeros(size(beta,1), nsim, length(h1));
mspe_mat_hn = zeros(size(beta,1), nsim, length(h2));
msee_mat_hn = zeros(size(beta,1), nsim, length(h2));

% initializing sum over nsim simulations~
mspe_matsum_ln = real(mspe_ln{1});
msee_matsum_ln = real(msee_ln{1});
mspe_matsum_hn = real(mspe_hn{1});
msee_matsum_hn = real(msee_hn{1});

for i=1:nsim
    for r=1:length(beta)
        for l=1:length(h1)
            mspe_mat_ln(r, i, l) = real(mspe_ln{i}(r,l));
            msee_mat_ln(r, i, l) = real(msee_ln{i}(r,l));
        end
        for l=1:length(h2)
            mspe_mat_hn(r, i, l) = real(mspe_hn{i}(r,l));
            msee_mat_hn(r, i, l) = real(msee_hn{i}(r,l));
        end
    end
end

% finding out sum and average of MSPE, MSEE over simulations for given theta, bandwidth  
for i=2:nsim
    
    mspe_matsum_hn = mspe_matsum_hn + real(mspe_hn{i});
    msee_matsum_hn = msee_matsum_hn + real(msee_hn{i});
    mspe_matsum_ln = mspe_matsum_ln + real(mspe_ln{i});
    msee_matsum_ln = msee_matsum_ln + real(msee_ln{i});
    
end

% averaging the MSPE/MSEE values over all nsim simulations
mspe_matsum_hn = mspe_matsum_hn./nsim;
msee_matsum_hn = msee_matsum_hn./nsim;
mspe_matsum_ln = mspe_matsum_ln./nsim;
msee_matsum_ln = msee_matsum_ln./nsim;

[x,y]=find(mspe_matsum_hn == min(mspe_matsum_hn(:)));
mspe_matsum_hn(x,y)
betaopt_mspe_hn= beta(x,:); % optimum theta that minimizes MSPE
h_opt_mspe_hn = h2(y); % optimum bandwidth that minimizes MSPE

[x,y]=find(msee_matsum_hn == min(msee_matsum_hn(:)));
msee_matsum_hn(x,y)
betaopt_msee_hn=  beta(x,:); % optimum theta that minimizes MSEE
h_opt_msee_hn= h2(y); % optimum bandwidth that minimizes MSEE

[x,y]=find(mspe_matsum_ln == min(mspe_matsum_ln(:)));
mspe_matsum_ln(x,y)
betaopt_mspe_ln= beta(x,:); % optimum theta that minimizes MSPE
h_opt_mspe_ln= h1(y); % optimum bandwidth that minimizes MSPE

[x,y]=find(msee_matsum_ln == min(msee_matsum_ln(:)));
msee_matsum_ln(x,y)
betaopt_msee_ln = beta(x,:); % optimum theta that minimizes MSEE
h_opt_msee_ln = h1(y); % optimum bandwidth that minimizes MSEE

distv= zeros(nsim, p);
rw1 = find(h1==h_opt_msee_ln);
rw2 = find(h2==h_opt_msee_hn);

for i=1:nsim
    distv(i, :)= beta(find(msee_mat_hn(:,i,rw2)==min(msee_mat_hn(:,i,rw2))),:);
end

distvs = distv(:,1:2);

















%% Sphere simulations n=50, p=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Adding path to files and folder
addpath(genpath('manopt'))

nsim = 5;  % number of simulations
n=200;        % number of observations 
p=2;         % dimension of theta vector 
s=19878;     % random number generators
rng(s);      % setting seed for random number generator

%Variances for normal data simulation
tau1=0.2;   % Low noise 
tau2=0.6;  % High Noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true theta
b =[cos(-0.704), sin(-0.704)];

% Creating the grid for theta on hypersphere p=4, using polar coordinates
f = -1.57:(pi/60):1.57;
beta = [cos(f)', sin(f)'];
beta = beta(beta(:,1)>0.02 ,:);  
beta = beta(beta(:,1)<0.98 ,:);  

h1 = exp(log(0.5):0.02:log(.85));
h2 = exp(log(0.65):0.02:log(1.5));

% h_opt_mspe_hn = .8438;   
% h_opt_mspe_ln = 0.5327;  
% h_opt_msee_hn = 0.8271;
% h_opt_msee_ln = .5434;

pm = struct('s', s, 'nsim', nsim, 'n', n, 'tau1', tau1,'tau2', tau2, 'h1', h1, 'h2', h2);

% Generating the data from fixed theta value, here we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

% Generate N(0, (0.2)^2) random numbers in matrix n x nsim in 2 sets u1 and u2 respectively
  u1 = tau1*randn(n, nsim); 
  u2 = tau1*randn(n, nsim); 
  
  % Generate N(0, (0.6)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  v1 = tau2*randn(n, nsim); 
  v2 = tau2*randn(n, nsim); 
  
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
mspe_ln = arrayfun(@(i) zeros(length(beta), length(h1)), 1:nsim, 'UniformOutput', false);
mspe_hn = arrayfun(@(i) zeros(length(beta), length(h2)), 1:nsim, 'UniformOutput', false);

% Creating the array with necessary dimensions to store MSEE in both high and low noises respectively
msee_ln = arrayfun(@(i) zeros(length(beta), length(h1)), 1:nsim, 'UniformOutput', false);
msee_hn = arrayfun(@(i) zeros(length(beta), length(h2)), 1:nsim, 'UniformOutput', false);

pool=parpool(6);  % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim % starting the parallel for loop
%for i = 1:nsim
disp(i)

for r = 1:length(beta) % starting for loop on the theta parameter.
    
disp(r)
xtemp = [x(:,i,1), x(:,i,2)];     % getting the covariates for the i-th simulation
ztemp = (xtemp*beta(r,:)')./sqrt(p);  % dividing by sqrt(p) so that the entries in ztemp in [-1,1].

    for l = 1:length(h1)
        
       Y_fr_lo = get_sphere_fit_SIM2(Y_true1{i}, ztemp, h1(l)); % getting the frechet mean for given theta, covariates
       mspe_ln{i}(r,l) = mean((arrayfun(@(g) acos(Y_true1{i}(:, g)'*Y_fr_lo(:,g)),1:n)).^2); % compute the MSPE
       msee_ln{i}(r,l) = mean((arrayfun(@(g) acos(mreg1{i}(:, g)'*Y_fr_lo(:,g)),1:n)).^2);    % compute the MSEE
    
    end
    for i2 = 1:length(h2)

       Y_fr_hi = get_sphere_fit_SIM2(Y_true2{i}, ztemp, h2(i2)); % getting the frechet mean for given theta, covariates
       mspe_hn{i}(r,i2) = mean((arrayfun(@(g) acos(Y_true2{i}(:, g)'*Y_fr_hi(:,g)),1:n)).^2); % compute the MSPE
       msee_hn{i}(r,i2) = mean((arrayfun(@(g) acos(mreg1{i}(:, g)'*Y_fr_hi(:,g)),1:n)).^2);    % compute the MSEE

    end
  
end

end
delete(pool) % shutting down parallel pool 

save Sphere_results_n200_nsim5_p2.mat mreg1 * mspe_hn * msee_hn * mspe_ln * msee_ln * h1 * h2

%the following are some computation necessary for further analyses~

mspe_mat_ln = zeros(size(beta,1), nsim, length(h1));
msee_mat_ln = zeros(size(beta,1), nsim, length(h1));
mspe_mat_hn = zeros(size(beta,1), nsim, length(h2));
msee_mat_hn = zeros(size(beta,1), nsim, length(h2));

% initializing sum over nsim simulations~
mspe_matsum_ln = real(mspe_ln{1});
msee_matsum_ln = real(msee_ln{1});
mspe_matsum_hn = real(mspe_hn{1});
msee_matsum_hn = real(msee_hn{1});

for i=1:nsim
    for r=1:length(beta)
        for l=1:length(h1)
            mspe_mat_ln(r, i, l) = real(mspe_ln{i}(r,l));
            msee_mat_ln(r, i, l) = real(msee_ln{i}(r,l));
        end
        for l=1:length(h2)
            mspe_mat_hn(r, i, l) = real(mspe_hn{i}(r,l));
            msee_mat_hn(r, i, l) = real(msee_hn{i}(r,l));
        end
    end
end

% finding out sum and average of MSPE, MSEE over simulations for given theta, bandwidth  
for i=2:nsim
    
    mspe_matsum_hn = mspe_matsum_hn + real(mspe_hn{i});
    msee_matsum_hn = msee_matsum_hn + real(msee_hn{i});
    mspe_matsum_ln = mspe_matsum_ln + real(mspe_ln{i});
    msee_matsum_ln = msee_matsum_ln + real(msee_ln{i});
    
end

% averaging the MSPE/MSEE values over all nsim simulations
mspe_matsum_hn = mspe_matsum_hn./nsim;
msee_matsum_hn = msee_matsum_hn./nsim;
mspe_matsum_ln = mspe_matsum_ln./nsim;
msee_matsum_ln = msee_matsum_ln./nsim;


[x,y]=find(mspe_matsum_hn == min(mspe_matsum_hn(:)));
mspe_matsum_hn(x,y)
betaopt_mspe_hn= beta(x,:); % optimum theta that minimizes MSPE
h_opt_mspe_hn = h2(y); % optimum bandwidth that minimizes MSPE

[x,y]=find(msee_matsum_hn == min(msee_matsum_hn(:)));
msee_matsum_hn(x,y)
betaopt_msee_hn=  beta(x,:); % optimum theta that minimizes MSEE
h_opt_msee_hn= h2(y); % optimum bandwidth that minimizes MSEE

[x,y]=find(mspe_matsum_ln == min(mspe_matsum_ln(:)));
mspe_matsum_ln(x,y)
betaopt_mspe_ln= beta(x,:); % optimum theta that minimizes MSPE
h_opt_mspe_ln= h1(y); % optimum bandwidth that minimizes MSPE

[x,y]=find(msee_matsum_ln == min(msee_matsum_ln(:)));
msee_matsum_ln(x,y)
betaopt_msee_ln = beta(x,:); % optimum theta that minimizes MSEE
h_opt_msee_ln = h1(y); % optimum bandwidth that minimizes MSEE

distv= zeros(nsim, p);
rw1 = find(h1==h_opt_msee_ln);
rw2 = find(h2==h_opt_msee_hn);

for i=1:nsim
    distv(i, :)= beta(find(msee_mat_hn(:,i,rw2)==min(msee_mat_hn(:,i,rw2))),:);
end

distvs = distv(:,1:2);
