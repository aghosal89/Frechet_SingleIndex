%% Computation of MSPE/MSEE for various simulation scenarios using Nelder-Mead optimization

%% Sphere simulations various n, p=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
% Adding path to 'manopt' folder
addpath(genpath('manopt'))
addpath(genpath('FMINSEARCHBND'))

% For beta in circumference of a unit circle;

nsim=200;    % number of simulations
n=200;        % number of observations 
p=2;         % dimension of theta vector 
s=19878;     % seed for n=50
%s=19879;     % seed for n=100
%s=19880;     % seed for n=200
rng(s);      % setting seed for random number generator

% Variances for normal data simulation
tau=0.6;    % Noise

%function to project a point in [-1,1] onto circumference of unit sphere with center at origin
reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

% Choose the true theta
if(p==2)
    b= [cos(-.704) sin(.704)];
elseif(p==3)
    b=[0.386, 0.181, -sqrt(1-0.386^2-0.181^2)];
elseif(p==4)
    b =[0.6087, 0.1332, -0.287, -sqrt(1-0.6087^2-0.1332^2-0.287^2)];
end

% Creating the starting points for eta on hypersphere p=2, store in X
nsp=5; % number of starting points 
spc= pi/nsp;  
f = linspace(-pi/2+spc/2,pi/2-spc/2,((pi-spc)/spc)+1); %choosing the increments

if(p==2)
    [S1] = ndgrid(f'); % creating the grid
    X = [S1(:)];
elseif(p==3)
    [S1 S2] = ndgrid(f', f'); % creating the grid
    X = [S1(:) S2(:)];
elseif(p==4)
    [S1 S2 S3] = ndgrid(f', f', f'); % creating the grid
    X = [S1(:) S2(:) S3(:)];
end

start=zeros(size(X,1), size(X,2)+1);
for i=1:size(start,1)
    start(i,:)= polar2cart(X(i,:),1);
end    

% Generating the data (x,Y) from fixed \theta_0 value, here we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 
  
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

  % Map the point in (-1,1) on to the circumference of unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false); 
  
  % spherefactory is a function in the 'manopt' that returns a manifold struct to optimize over unit-norm vectors or matrices.
  % defining a set of functions pre-defined on a unit sphere of p-dimensions
  M = spherefactory(3); % define structures fr, ops to store the pre-defined functions on a sphere
  fr = struct(); 
  ops = struct();
  fr.M = M; 
  ops.verbosity = 0;
 
  % create the variables for storing the simulated response data based on
  % our assumptions on theta and covariates: creates a list of nsim number
  % of matrices of size 3 x n, we have 
  
  Y = arrayfun(@(l) zeros(3, size(z,1)), 1:nsim, 'UniformOutput', false);
  
  % for each simulation, generate response data by first projecting the
  % single index z (a nx1 matrix)
  % onto the unit sphere, then map those n points to the tangent space, and generate a single point 
  % by adding noise created by normal samples in v1, v2.
  
  for i=1:nsim
     Y{i} = cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j),[v1(j, i) v2(j, i)]), 1:size(z,1),'UniformOutput', false)); 
  end
  
% choosing a range bandwidths~
bw_min= zeros(size(start,1), nsim);
bw_max= zeros(size(start,1), nsim);
temp=zeros(nsim,1);
for k=1:size(start,1)
    for i=1:nsim
        
        if(p==2)
            temp =(start(k,1)*x(:,i,1) + start(k,2)*x(:,i,2))./sqrt(p);
        elseif(p==3)
            temp =(start(k,1)*x(:,i,1) + start(k,2)*x(:,i,2)+ start(k,3)*x(:,i,3))./sqrt(p);
        elseif(p==4)
            temp =(start(k,1)*x(:,i,1) + start(k,2)*x(:,i,2)+start(k,3)*x(:,i,3)+start(k,4)*x(:,i,4))./sqrt(p);
        end
        
        tempd=abs(repmat(temp, 1, length(temp)) - repmat(temp', length(temp), 1));
        tempc=sort(tempd,2);
        bw_min(k,i)= max(tempc(:,4));
        bw_max(k,i)=max(tempd(:));
    end
end    

h=exp(linspace(log(max(bw_min(:))*1.25),log(max(bw_max(:))),10));

pm = struct('s', s, 'nsim', nsim, 'n', n, 'h', h, 'tau', tau);

% Creating the array with necessary dimensions to store MSPE, MSEE in both high and low noises respectively 

mspe = arrayfun(@(i) zeros(length(h), size(start,1)), 1:nsim, 'UniformOutput', false);
msee = arrayfun(@(i) zeros(length(h), size(start,1)), 1:nsim, 'UniformOutput', false);

eta_opt = arrayfun(@(l) zeros(length(h), p-1, size(start,1)), 1:nsim, 'UniformOutput', false);

%pool=parpool(6);   % starting computation in a parallel pool with 6 workers
%parfor i = 1:nsim  % starting the parallel for loop
for i = 1:nsim
disp(i)

Y_temp = Y{i};             % creating temporary response variable to be input for cost function

% getting the covariates for the i-th dataset
if(p==2)
    xtemp = [x(:,i,1) x(:,i,2)];
elseif(p==3)
    xtemp = [x(:,i,1) x(:,i,2) x(:,i,3)];
elseif(p==4)
    xtemp = [x(:,i,1) x(:,i,2) x(:,i,3) x(:,i,4)];
end 
              
for j=1:size(start,1)

    for l = 1:length(h) % bandwidth range for high noise scenario
    disp(['NM_nsim200_n200_p2_HN' '_' 'nsim' '_' num2str(i, '%03d') '_' 'start' '_' num2str(j, '%03d') '_' 'bandwidth' '_' num2str(l, '%02d')])
    
    cost_fun = @(eta) nm_cost(Y_temp, eta, h(l), xtemp); % define cost function 
       
    % running the optimization for starting point 
    
    if(p==2)
        [eta_opt1, fval] = fminsearchbnd(cost_fun, X(j,:), -pi/2, pi/2);
    elseif(p==3)
        [eta_opt1, fval] = fminsearchbnd(cost_fun, X(j,:), [-pi/2,-pi/2],[pi/2, pi/2]);
    elseif(p==4)
        [eta_opt1, fval] = fminsearchbnd(cost_fun, X(j,:), [-pi/2,-pi/2,-pi/2],[pi/2, pi/2,pi/2]);
    end
    
    mspe{i}(l, j) = fval;           % storing minimum cost function
    
    eta_opt{i}(l,:,j) = eta_opt_temp;   % storing local minima
  
    end
  
end

end

%delete(pool) % shutting down parallel pool 

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
mspe_nm_sum = sum(mspe_nm,2);
h_opt=h(find(mspe_nm_sum == min(mspe_nm_sum(:))));

% mean and standard deviation of the MSPE over simulations~
mean(mspe_nm(find(mspe_nm_sum == min(mspe_nm_sum(:))), :));
std(mspe_nm(find(mspe_nm_sum == min(mspe_nm_sum(:))), :));

% distribution of beta estimates 
eta_opt_s = zeros(nsim, p-1);
beta_opt_s = zeros(nsim, p);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
    beta_opt_s(i,:)= polar2cart(eta_opt_s(i,:),1);
end


% MSE  of theta estimator
mean(acos(beta_opt_s*b'))
std(acos(beta_opt_s*b'))
 
% computing MSEE from the estimated beta coefficients~
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

save NM_Sphere_results_n200_nsim200_p2_HN.mat * 

