%% This script is for testing out the quality of the gradient and Hessian approximations in WnCost.m and getWnGradHess.m

clear all; close all;
% Adding path to access files in 'manopt' folder
addpath(genpath('manopt'))
% Adding path to access files in 'FMINSEARCHBND' folder for Nelder-Mead optimization
addpath(genpath('FMINSEARCHBND'))


%% Data Generation Settings
n = 20;        % number of observations 
p = 2;         % dimension of index parameter 
noise = 'Low'; % High or low noise setting

% Parameters used to initialize the optimization for theta
initType = 'fixed'; % 'random' or 'fixed' - how should the initialization grid be chosen?
nsp = 3; % number of starting points to be used in each angular component, if initType == 'fixed'
nrnd = 3; 

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
theta0 = tmp/norm(tmp); % normalize to have norm 1 (this is theta_0 in the paper)

%% Generate Data

rng(s_dt);      % setting seed for generating data
  
% Generate Uniform(0,1) random variables in matrix n x nsim x p
  % dimensions, they are used as covariates, p gives us the number of
  % covariates used in the model. Hence the data space is (0,1)^p
  
  x = rand(n, p);  % generate the explanatory variables
  
% generate the index (linear combinations) with given parameters, this
  % gives us the matrix of dimension n x nsim , we divide by sqrt(p) so that all elements of z 
  % are in the range (-1,1), as required by the function mreg below.  
  
  z = x*theta0/sqrt(p);
  
% Generating the data (x,Y) from fixed index parameter value above for a given p, we assume that the
% parameter space is a proper subset of a p-dimensional unit sphere with
% center at origin and first component positive. 

  % Map the point in (-1,1) on to the circumference of unit sphere, this
  % gives us the conditional expectation of Y given x
  mreg = reg_curve(z)'; 

  % Generate N(0, (tau)^2) random numbers in matrix n x nsim in 2 sets v1 and v2 respectively
  % These represent random errors in the tangent space 
  v1 = tau*randn(n, 1); 
  v2 = tau*randn(n, 1); 

 
  % create the object for storing the simulated response data based on
  % our assumptions on theta and covariates: 
  
  Y = cell2mat(arrayfun(@(j) add_noise(mreg(:, j), [v1(j) v2(j)]), 1:n, 'UniformOutput', false));
  
%% Create grid of starting values for theta and grids of bandwidths

% Creating the grid of polar coordinates as starting points for
% p-dimensional hypersphere, stored in theta_init

if(strcmp(initType, 'fixed'))
    theta_init = thetaGrid(nsp, p - 1);
else    
    theta_init = nrnd;
end

% choose a bandwidth

bw_min = 0;
bw_max = 0;

    normArr = zeros(n, n);
    for i = 1:(n - 1)
        for j = (i + 1):n
            normArr(i, j) = norm(x(i, :) - x(j, :));
            normArr(j, i) = normArr(i, j);
        end
    end

    normSt = sort(normArr, 2);
    bw_min = max(bw_min, max(normSt(:, 4))); bw_max = max(bw_max, max(normSt(:, n)));

% make bandwidth as geometric mean of min and max
h = exp((log(bw_min) + log(bw_max))/2);

%% This chunk is for comparing Nelder-Mead to the new trust regions optimizer - only run if derivatives are working fine
clear optns
optns(6) = struct();

optns(1).Algorithm = 'NM'; optns(1).Ret = 0; % Should run Nelder-Mead with 100 max iterations
optns(2).Algorithm = 'NM'; optns(2).Ret = 0; optns(2).useAltCost = 1; % NM with alternative cost function
optns(3).Algorithm = 'TR'; optns(3).Ret = 0; optns(3).useAltGrad = 1; % TR with gradient approximation
optns(4).Algorithm = 'TR'; optns(4).Ret = 0; optns(4).useAltGrad = 1; optns(4).useAltHess = 1; % TR with gradient and Hessian approximation
optns(5).Algorithm = 'TR'; optns(5).Ret = 0; optns(5).useAltCost = 1; optns(5).useAltGrad = 1; % TR with cost and gradient approximation
optns(6).Algorithm = 'TR'; optns(6).Ret = 0; optns(6).useAltCost = 1; optns(6).useAltGrad = 1; optns(6).useAltHess = 1; % TR with cost, gradient and Hessian approximation

fsi_fits = cell(1, 6);
tms = zeros(1, 6);

for j = 1:6
    
    disp(j)
    tic;
    fsi_fit{j} = get_sphere_fit_FSI(Y, x, h, [], theta_init, optns(j));
    tms(j) = toc;
    
end

%% Evaluate cost function for all theta values in the grid and also the gradient and Hessian approximations

Wn = zeros(1, nsp);
WnGrad = zeros(1, nsp);
WnHess = zeros(1, nsp);
WnCostAlt = zeros(1, nsp);
eta = arrayfun(@(j) cart2polar(theta_init(:, j)), 1:nsp);

for j = 1:nsp
    disp(j)
    [Wn(j), WnGrad(j), WnHess(j), WnCostAlt(j)] = WnCost(eta(j), Y, x, h);
    
end


%% Some derivative checks
eta = arrayfun(@(j) cart2polar(theta_init(:, j)), 1:size(theta_init, 2));
del = eta(2) - eta(1);
mu0Mat = zeros(n, size(theta_init, 2)); mu1Mat = mu0Mat; mu2Mat = mu0Mat; sig0sqMat = mu0Mat;
h0Mat = mu0Mat; h1Mat = mu0Mat; h2Mat = mu0Mat; lMat = mu0Mat; % compare these to the numeric derivatives of mu0Mat, mu1Mat, mu2Mat, and sig0sqMat
H0Mat = mu0Mat; H1Mat = mu0Mat; h2Mat = mu0Mat; LMat = mu0Mat; % hessians of the above

vArr = zeros(n, n, size(theta_init, 2)); mArr = vArr; wArr = vArr; qArr = vArr; %  mArr is the deriv of vArr and qArr is the deriv of wArr
MArr = vArr; QArr = vArr; % Hessians

YstarArr = zeros(3, n, size(theta_init, 2)); YtilArr = YstarArr; % Local linear fits and normalized LL fits 
tArr = zeros(n, p - 1, 3, size(theta_init, 2)); ttilArr = tArr; % Derivs of local linear fits and normalized LL fits
TArr = tArr; TtilArr = tArr; % Hessians

WnVec = zeros(1, size(theta_init, 2)); WnDerVec = WnVec; WnHessVec = WnVec; % approximate criterion values and their derivatives

for s = 1:size(theta_init, 2)
    
    J = zeros(p - 1, p); % Jacobian matrix
    etaA = [eta(s); pi/2]; % augmented eta vector

    for j = 1:(p - 1)
        for k = 1:p
        
            if j > (p - k + 1)    
                J(j, k) = 0; % eta_j not involved in computing theta_k
            elseif j == (p - k + 1)
                J(j, k) = prod(cos(etaA(1:(p - k + 1))));
            else
                tmp = cos(etaA); tmp(j) = sin(etaA(j)); tmp(p - k + 1) = sin(etaA(p - k + 1));
                J(j, k) = -prod(tmp);
            end
        end
    end
    
    grad = 0; % Incrementally compute the gradient
    hess = 0;
    
    for i = 1:n
    
        zcur = x*theta_init(:, s);
        zdiff = zcur - zcur(i);    
        xdiff = x - repmat(x(i, :), n, 1);
    
        % Get normalized local linear fits (these hopefully approximate the LF
        % fits)
 
            wcur = getLFRweights(zcur, zcur(i), h);
            Ystar = sum(cell2mat(arrayfun(@(j) wcur(j)*Y(:, j), 1:n, 'UniformOutput', false)), 2);
            Ytil = Ystar/norm(Ystar);
            wArr(:, i, s) = wcur; YstarArr(:, i, s) = Ystar; YtilArr(:, i, s) = Ytil; % store values
            
        % Get derivatives of local linear weights at x(i, :)
        
        % Kernel weights/derivatives and averages associated with local
        % linear estimation
        
            [Kvec, KvecDer, KvecDer2] = K(zdiff, h);
            mu0 = mean(Kvec);  mu1 = mean(Kvec.*zdiff); mu2 = mean(Kvec.*(zdiff.^2));
            sig0sq = mu0*mu2 - mu1^2;
            mu0Mat(i, s) = mu0; mu1Mat(i, s) = mu1; mu2Mat(i, s) = mu2; sig0sqMat(i, s) = sig0sq; % store values
    
        % Derivatives (w.r.t. theta) of muj (for fitting at x(i, :))
        
            h0 = mean(repmat(KvecDer, 1, p).*xdiff)';
            h1 = mean(repmat(KvecDer.*zdiff + Kvec, 1, p).*xdiff)'; 
            h2 = mean(repmat(KvecDer.*(zdiff.^2) + 2*Kvec.*zdiff, 1, p).*xdiff)';
            h0Mat(i, s) = J*h0; h1Mat(i, s) = J*h1; h2Mat(i, s) = J*h2; % store derivs w.r.t. eta
            
        % More derivatives
        
            l = mu2*h0 + mu0*h2 - 2*mu1*h1; % deriv of sig0sq
            v = mu2 - mu1*zdiff; 
            m = repmat(h2, 1, n) - h1*zdiff' - mu1*xdiff'; % derivs of v = (mu2 - m1*zdiff)
            u = repmat((KvecDer.*v)', p, 1).*xdiff' + repmat(Kvec', p, 1).*m; % derivs of numerator of local linear weight
            q = (u/sig0sq - l*(Kvec.*v)'/(sig0sq^2))/n; % derivs of local linear weights
            lMat(i, s) = J*l; vArr(:, i, s) = v; mArr(:, i, s) = m'*J'; qArr(:, i, s) = q'*J'; % store values
            
            
        % Get derivative of normalized local linear fit at x(i, :)    
        
            t = 0;
            for j = 1:n
                t = t + q(:, j)*Y(:, j)';
            end % t = derivative of local linear fit at x(i, :)
    
            b = (eye(3) - Ystar*Ystar'/(norm(Ystar)^2))/norm(Ystar); % derivative of normalization
            ttil = t*b; % derivative of normalized local linear fit at x(i, :)
            tArr(i, :, :, s) = J*t; ttilArr(i, :, :, s) = J*ttil; % store values
            
            ip = Y(:, i)'*Ytil;
            grad = grad - 2*acos(ip)*(1 - ip^2)^(-1/2)*ttil*Y(:, i);
            
            H0 = zeros(p, p); H1 = H0; H2 = H0;
                
            for j = 1:n
            
                H0 = H0 + KvecDer2(j)*xdiff(j, :)'*xdiff(j, :)/n;
                H1 = H1 + (KvecDer2(j)*zdiff(j) + 2*KvecDer(j))*xdiff(j, :)'*xdiff(j, :)/n;
                H2 = H2 + (KvecDer2(j)*zdiff(j)^2 + 4*KvecDer(j)*zdiff(j) + 2*Kvec(j))*xdiff(j, :)'*xdiff(j, :)/n;
                
            end
            H0Mat(i, s) = J*H0*J'; H1Mat(i, s) = J*H1*J'; H2Mat(i, s) = J*H2*J'; % store values
            
        % More Hessians
        
            L = h0*h2' + mu2*H0 + h2*h0' + mu0*H2 - 2*h1*h1' - 2*mu1*H1; % Hess of sig0sq
            M = zeros(p, p, n); Q = M;
            for j = 1:n
                M(:, :, j) = H2 - h1*xdiff(j, :) - zdiff(j)*H1 - xdiff(j, :)'*(h1'); % store values
                Uj = KvecDer2(j)*v(j)*xdiff(j, :)'*xdiff(j, :) + KvecDer(j)*(xdiff(j, :)'*m(:, j)' + m(:, j)*xdiff(j, :)) + Kvec(j)*squeeze(M(:, :, j));
                Q(:, :, j) = (Uj/sig0sq - u(:, j)*l'/(sig0sq^2) - l*u(:, j)'/(sig0sq^2) - Kvec(j)*v(j)*L/(sig0sq^2) + 2*Kvec(j)*v(j)*l*l'/(sig0sq^3))/n; % store values
            end % Hessians of v = (mu2 - m1*zdiff)
            
            LMat(i, s) = J*L*J'; MArr(:, i, s) = arrayfun(@(j) J*squeeze(M(:, :, j))*J', 1:n); QArr(:, i, s) = arrayfun(@(j) J*squeeze(Q(:, :, j))*J', 1:n); % store values - the second two commands will only work when p = 2, need to be generalized for larger p
            
        % Get Hessian of normalized local linear fit at x(i, :)
    
            T = zeros(p, 3, p);
            for j = 1:n
                for k = 1:p
                    T(:, :, k) = T(:, :, k) + squeeze(Q(:, k, j))*Y(:, j)';
                end
            end % Hessian of local linear fit at x(:, i)
            TArr(i, :, :, s) = arrayfun(@(l) J*squeeze(T(:, l, :))*J', 1:3); % store values
            
            B = zeros(3, 3, 3);
            for k = 1:3
                ek = zeros(3, 1); ek(k) = 1;
                B(:, :, k) = -Ystar(k)*eye(3)/(norm(Ystar)^3) - (ek*Ystar' + Ystar*ek')/(norm(Ystar)^3) + 3*Ystar*Ystar'*Ystar(k)/(norm(Ystar)^5);
            end % Hessian of normalization
        
            Ttil = zeros(p, 3, p);
            for k = 1:p
            
                tmp = zeros(3, 3);
                for kk = 1:3
                    for ll = 1:3
                        tmp(kk, ll) = t(k, :)*squeeze(B(kk, ll, :));
                    end
                end
            
                Ttil(:, :, k) = squeeze(T(:, :, k))*b + t*tmp;
            
            end
            % Hessian of normalized local linear fit at x(:, i)
            TtilArr(i, :, :, s) = arrayfun(@(l) J*squeeze(Ttil(:, l, :))*J', 1:3); % store values
            
            tmp = cell2mat(arrayfun(@(k) squeeze(Ttil(:, :, k)*Y(:, i)), 1:p, 'UniformOutput', false));
            hess = hess - 2*(ttil*Y(:, i)*Y(:, i)'*ttil')*(1 - ip^2)^(-1)*(1 + acos(ip)*(1 - ip^2)^(-1/2)) + acos(ip)*(1 - ip^2)^(-1/2)*tmp;
            hess = (hess + hess')/2;
            
    end
    
    WnVec(s) =  sum((arrayfun(@(g) acos(Y(:, g)'*YtilArr(:, g, s)), 1:n)).^2); % prediction error on sphere
    WnDerVec(s) = J*grad;
    WnHessVec(s) = J*hess*J';
    
end

% Plots
% Fix x and z in formula
i = 10;
ii = 15;

% h plots
figure; hold on;
plot(eta, h0Mat(i, :), 'k-')
plot(eta, gradient(mu0Mat(i, :), del), 'ro')

figure; hold on;
plot(eta, h1Mat(i, :), 'k-')
plot(eta, gradient(mu1Mat(i, :), del), 'ro')

figure; hold on;
plot(eta, h2Mat(i, :), 'k-')
plot(eta, gradient(mu2Mat(i, :), del), 'ro')

% m and q plots
figure; hold on;
plot(eta, squeeze(mArr(ii, i, :)), 'k-')
plot(eta, gradient(squeeze(vArr(ii, i, :)), del), 'ro')

figure; hold on;
plot(eta, squeeze(qArr(ii, i, :)), 'k-')
plot(eta, gradient(squeeze(wArr(ii, i, :)), del), 'ro')

% t plots

for j = 1:3
    
    figure; hold on;
    plot(eta, squeeze(tArr(i, :, j, :)), 'k-')
    plot(eta, gradient(squeeze(YstarArr(j, i, :)), del), 'ro')

end

% ttil plots

for j = 1:3
    
    figure; hold on;
    plot(eta, squeeze(ttilArr(i, :, j, :)), 'k-')
    plot(eta, gradient(squeeze(YtilArr(j, i, :)), del), 'ro')

end


% H plots

figure; hold on;
plot(eta, H0Mat(i, :), 'k-')
plot(eta, gradient(h0Mat(i, :), del), 'ro')

figure; hold on;
plot(eta, H1Mat(i, :), 'k-')
plot(eta, gradient(h1Mat(i, :), del), 'ro')

figure; hold on;
plot(eta, H2Mat(i, :), 'k-')
plot(eta, gradient(h2Mat(i, :), del), 'ro')


% M and Q plots
figure; hold on;
plot(eta, squeeze(MArr(ii, i, :)), 'k-')
plot(eta, gradient(squeeze(mArr(ii, i, :)), del), 'ro')

figure; hold on;
plot(eta, squeeze(QArr(ii, i, :)), 'k-')
plot(eta, gradient(squeeze(qArr(ii, i, :)), del), 'ro')

% T plots

for j = 1:3
    
    figure; hold on;
    plot(eta, squeeze(TArr(i, :, j, :)), 'k-')
    plot(eta, gradient(squeeze(tArr(i, :, j, :)), del), 'ro')

end

% Ttil plots

for j = 1:3
    
    figure; hold on;
    plot(eta, squeeze(TtilArr(i, :, j, :)), 'k-')
    plot(eta, gradient(squeeze(ttilArr(i, :, j, :)), del), 'ro')

end


% Wn plots

figure; hold on
plot(eta, WnDerVec, 'k-')
plot(eta, gradient(WnVec, del), 'ro')

figure; hold on;
plot(eta, WnHessVec, 'k-')
plot(eta, gradient(WnDerVec, del), 'ro')