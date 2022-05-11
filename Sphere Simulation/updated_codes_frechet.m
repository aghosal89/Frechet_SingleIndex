%% Computation of MSPE/MSEE for various simulation scenarios using full grid-search optimization

%% Sphere simulations for various n, p 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;   

% Adding path to access files in 'manopt' folder
addpath(genpath('manopt'))

nsim = 200;   % number of simulations
n = 200;      % number of observations 
p = 4;        % dimension of index parameter 
s=19378;      % seed 
rng(s);       % setting seed for random number generator

%Variances for data simulation
tau=0.6;  % Noise level

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true theta for the number of covariates p:
if(p==2)
    b = [cos(-.704) sin(.704)];
elseif(p==3)
    b = [0.386, 0.181, -sqrt(1-0.386^2-0.181^2)];
elseif(p==4)
    b = [0.7430, 0.0120, -0.4870, 0.4590];
end

% Creating the grid of polar coordinates as the full grid points for p-dimensional hypersphere, 

% number of starting points in each angular component as varied with p
if(p==2)
    nsp=100; 
elseif(p==3)
    nsp=20;
elseif(p==4)
    nsp=10;
end

spc= pi/nsp;  % increment
f = linspace(-pi/2+spc/2,pi/2-spc/2,((pi-spc)/spc)+1); % equidistant starting points

% the full grid stored stored in X
if(p==2)
    [S1] = ndgrid(f'); % creating the grid for p=2
    X = [S1(:)];
elseif(p==3)
    [S1 S2] = ndgrid(f', f'); % creating the grid for p=3
    X = [S1(:) S2(:)];
elseif(p==4)
    [S1 S2 S3] = ndgrid(f', f', f'); % creating the grid for p=4
    X = [S1(:) S2(:) S3(:)];
end

% converting polar coordinates of grid points into cartesian coordinates
beta=zeros(size(X,1), size(X,2)+1);
for i=1:size(beta,1)
    beta(i,:)= polar2cart(X(i,:),1);
end    

% Generate N(0, (tau)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  v1 = tau*randn(n, nsim); 
  v2 = tau*randn(n, nsim); 

% Generate Uniform(0,1) random variables in matrix n x nsim x p
  % dimensions, they are used as covariates, p gives us the number of
  % covariates used in the model. Hence the data space is (0,1)^p
  
  x = rand(n, nsim, p);  % generate the explanatory variables
% generate the index (linear combinations) with given parameters, this
  % gives us the matrix of dimension n x nsim , we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), now this z can be used as the explanatory variable in Local 
  % frechet regression, giving us the Frechet Single Index regression. 
  
  if(p==2)
      z= (b(1)*x(:,:,1) + b(2)*x(:,:,2))./sqrt(p); 
  elseif(p==3)
      z= (b(1)*x(:,:,1) + b(2)*x(:,:,2) + b(3)*x(:,:,3))./sqrt(p);
  elseif(p==4)
      z= (b(1)*x(:,:,1) + b(2)*x(:,:,2) + b(3)*x(:,:,3) + b(4)*x(:,:,4))./sqrt(p);
  end

% Generating the data (x,Y) from fixed index parameter value above for a given p, we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

  % Map the point in (-1,1) on to the circumference of unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
  % spherefactory is a function in the 'manopt' folder that returns a manifold struct to optimize over 
  % unit-norm vectors or matrices, defining a set of functions pre-defined on a unit sphere of p-dimensions
  
  M = spherefactory(3); % define structures fr, ops to store the pre-defined functions on a sphere
  fr = struct(); 
  ops = struct();
  fr.M = M; 
  ops.verbosity = 0;
 
  % create the object for storing the simulated response data based on
  % our assumptions on theta and covariates: 
  
  Y = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  
  % for each simulation, the response data is generated from mreg, the
  % projections onto the unit sphere, by mapping them to the tangent space,
  % and generate a single point (a 2-dimensional vector, sample from a bivariate normal distribution) 
  % by adding noise created by normal samples in v1, v2.
  % then we map it back to the surface of unit sphere and get the response.
  
  for i=1:nsim
     Y{i} = cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
% Creating the grid for theta on hypersphere, using polar coordinates
bw_min= zeros(size(beta,1),nsim);
bw_max= zeros(size(beta,1),nsim);

for k=1:size(beta,1)
    for i=1:nsim
        
        if(p==2)
            temp =(beta(k,1)*x(:,i,1)+beta(k,2)*x(:,i,2))./sqrt(p);
        elseif(p==3)
            temp =(beta(k,1)*x(:,i,1)+beta(k,2)*x(:,i,2)+ beta(k,3)*x(:,i,3))./sqrt(p);
        elseif(p==4)
            temp =(beta(k,1)*x(:,i,1)+beta(k,2)*x(:,i,2)+beta(k,3)*x(:,i,3)+beta(k,4)*x(:,i,4))./sqrt(p);
        end
        
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
% for i = 1:nsim  % starting the for loop
disp(i)

for r = 1:size(beta,1) % starting for loop on the grid of index parameters.
    
disp(r)

if(p==2)
    xtemp = [x(:,i,1) x(:,i,2)];
elseif(p==3)
    xtemp = [x(:,i,1) x(:,i,2) x(:,i,3)];
elseif(p==4)
    xtemp = [x(:,i,1) x(:,i,2) x(:,i,3) x(:,i,4)];
end 

% creating the temporary single index for our regression
% dividing single index by sqrt(p) so that the entries lie in [-1,1].
ztemp = (xtemp*beta(r,:)')./sqrt(p);  
  
    for l = 1:length(h)

       Y_fr = get_sphere_fit_FSI(Y{i}, ztemp, h(l)); % getting the frechet mean for given beta, Single index
       mspe{i}(r,l) = mean((arrayfun(@(g) acos(Y{i}(:, g)'*Y_fr(:,g)),1:n)).^2); % compute the MSPE
       msee{i}(r,l) = mean((arrayfun(@(g) acos(mreg{i}(:, g)'*Y_fr(:,g)),1:n)).^2); % compute the MSEE

    end
  
end

end
delete(pool) % shutting down parallel pool 

% The MSPE and MSEE values corresponding to each siulation, index
% parameter, bandwidth should be stored as matrix for further analysis.

mspe_mat = NaN(size(beta,1), nsim, length(h));
msee_mat = NaN(size(beta,1), nsim, length(h));

for i=1:nsim
    for r=1:length(beta)
        for q=1:length(h)
            mspe_mat(r, i, q) = (mspe{i}(r,q));
            msee_mat(r, i, q) = (mspe{i}(r,q));
        end
    end
end

% Computation for the average and standard deviation of the MSPE, MSEE over simulations
mspe_sum_mat = mspe{1}; % a matrix of dimension size(beta,1) x length(h)
msee_sum_mat = msee{1}; % a matrix of dimension size(beta,1) x length(h)

for i=2:nsim 
    mspe_sum_mat = mspe_sum_mat + mspe{i}; 
    msee_sum_mat = msee_sum_mat + msee{i}; 
end

mspe_ave_mat = mspe_sum_mat./nsim; % also a matrix of dimension size(beta,1) x length(h)

% computing the bandwidth that minimizes the MSPE
[~, c2] = find(mspe_ave_mat==min(mspe_ave_mat(:))); 
h_opt= h(c2);

% computing the estimate of single index for each dataset
beta_opt = zeros(nsim, size(beta,2));
for i=1:nsim
    beta_opt(i,:) = beta(find(mspe_mat(:,i,c2)==min(mspe_mat(:,i,c2))),:);
end    

% mean and standard error of the theta estimates
mean(acos(beta_opt*b')) 
std(acos(beta_opt*b'))

save Grid_Sphere_results_n200_nsim200_p4_HN.mat *




