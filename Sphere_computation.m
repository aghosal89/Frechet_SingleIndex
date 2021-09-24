
% For Nelder-Mead~
% For low noise
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
mspe_nm_sum = sum(mspe_nm,2)
h_opt=h(find(mspe_nm_sum == min(mspe_nm_sum(:))))

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
mean(acos(beta_opt_ar*b'))

% plot 2d histogram for p=3
addpath hist2d
hist2d(eta_opt_s(:,1), eta_opt_s(:,2), 'probability')

