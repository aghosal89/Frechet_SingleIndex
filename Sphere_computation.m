
% For Nelder-Mead~
% For low noise
mspe_nm = NaN(length(h), nsim);  % matrix to store the minimized MSPE for each bandwidth and each simulation
eta_hat = NaN(length(h), p-1, nsim);  % 3-D array to store the estimated index parameter corresponding to each bandwidth and each simulation.


for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); % for each bandwidth and each simulation find index providing which starting value provides the best estimate of index parameter
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2); % store estimate of eta corresponding to a bandwidth for each simulated data.
        mspe_nm(j,i) = mspe{i}(j, a2); % store minimized MSPE corresponding to a bandwidth for each simulated data.
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2) % sum of MSPE across all simulations 
h_opt=h(find(mspe_nm_sum == min(mspe_nm_sum(:)))) % find the best bandwidth

% mean and standard deviation of the MSPE over simulations~
mean(mspe_nm(find(mspe_nm_sum == min(mspe_nm_sum(:))), :)) 
std(mspe_nm(find(mspe_nm_sum == min(mspe_nm_sum(:))), :))

eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

% distribution of beta estimates 

beta_opt_s = zeros(nsim, p);
for i=1:nsim
    beta_opt_s(i,:)= polar2cart(eta_opt_s(i,:),1); 
end

% MSE  of theta estimator
mean(acos(beta_opt_ar*b')) % mean square error of theta estimates

% plot 2d histogram for p=3
addpath hist2d
hist2d(eta_opt_s(:,1), eta_opt_s(:,2), 'probability')

